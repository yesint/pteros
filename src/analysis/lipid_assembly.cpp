/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/


#include "pteros/analysis/lipid_assembly.h"
#include "pteros/core/pteros_error.h"
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include <Eigen/Dense>
#include <fstream>
#include <stack>
#include <sstream>

using namespace std;
using namespace pteros;
using namespace Eigen;

Lipid_assembly::Lipid_assembly(System &system, std::string lipids_sel, string heads_sel, string tails_sel, float d){
    create(system,lipids_sel,heads_sel,tails_sel,d);
}

string tcl_arrow(const Vector3f& p1, const Vector3f& p2, float r, string color){
    stringstream ss;
    Vector3f p = (p2-p1)*0.8+p1;
    ss << "draw color " << color << endl;

    ss << "draw cylinder \"" << p1.transpose()*10.0 << "\" ";
    ss << "\"" << p.transpose()*10.0 << "\" radius " << r << endl;

    ss << "draw cone \"" << p.transpose()*10.0 << "\" ";
    ss << "\"" << p2.transpose()*10.0 << "\" radius " << r*3.0 << endl;

    return ss.str();
}

void fit_quad_surface(const MatrixXf& coord,
                      Matrix<float,6,1>& res,
                      Vector3f& smoothed_first_point,
                      float& rms){
    int N = coord.cols();
    // We fit with polynomial fit = A*x^2 + B*y^2 + C*xy + D*x + E*y + F
    // Thus we need a linear system of size 6
    Matrix<float,6,6> m;
    Matrix<float,6,1> rhs; // Right hand side and result

    // Functions returning powers after coefficients (x^2, y^2, x*y, x, y, 1)
    vector<boost::function<float(int)> > coef(6);
    coef[0] = [&coord](int j){return coord.col(j)(0)*coord.col(j)(0);};
    coef[1] = [&coord](int j){return coord.col(j)(1)*coord.col(j)(1);};
    coef[2] = [&coord](int j){return coord.col(j)(0)*coord.col(j)(1);};
    coef[3] = [&coord](int j){return coord.col(j)(0);};
    coef[4] = [&coord](int j){return coord.col(j)(1);};
    coef[5] = [&coord](int j){return 1.0;};

    // Now form the matrix
    m.fill(0.0);
    rhs.fill(0.0);

    for(int r=0;r<6;++r){ //rows
        for(int c=0;c<6;++c){ //columns
            m(r,c) = 0.0;
            for(int j=0;j<N;++j){
                m(r,c) += coef[r](j)*coef[c](j);
            }
        }
        // Now rhs
        for(int j=0;j<N;++j){
            rhs(r) += coef[r](j)*coord.col(j)(2);
        }
    }

    // Now solve
    res = m.colPivHouseholderQr().solve(rhs);

    // Compute RMS of fitting
    rms = 0.0;
    for(int j=0;j<N;++j){
        float fit = 0.0;
        for(int r=0;r<6;++r) fit += coef[r](j)*res[r];
        rms += pow(coord.col(j)(2)-fit,2);
    }
    rms = sqrt(rms/float(N));

    // Get smoothed surface point
    smoothed_first_point.fill(0.0);
    for(int r=0;r<6;++r) smoothed_first_point(2) += coef[r](0)*res[r];
}

Local_properties Lipid_assembly::get_local_curvature(Selection& surf_spot, Selection& tail_spot){

    Local_properties prop;

    System* sys = surf_spot.get_system();
    int frame = surf_spot.get_frame();

    // Central marker atom (by convension the first atom in surf_spot)
    Vector3f marker_coord = surf_spot.XYZ(0);

    // Compute the inertia axes
    Vector3f moments;
    Matrix3f axes;
    surf_spot.inertia(moments,axes,true);
    // axes.col(2) is a local normal
    prop.surface_normal = axes.col(2).normalized();

    // Now compute the center of masses of tails ends
    Vector3f end = sys->Box(frame).get_closest_image(tail_spot.center(true,true),marker_coord);
    Vector3f lip_vec = marker_coord-end;
    // Now an angle between the lipid vector and normal
    float a = acos(prop.surface_normal.dot(lip_vec)/(prop.surface_normal.norm() * lip_vec.norm()));
    if(a>M_PI_2){ // Need to flip normal
        prop.surface_normal = -prop.surface_normal;
        axes.col(2) = -axes.col(2).eval();
    }

    // Now fit a quadric surface to local points

    // First compute transformation matrix to local basis of inertia axes
    Matrix3f tr,tr_inv;
    tr.col(0) = axes.col(0).normalized();
    tr.col(1) = axes.col(1).normalized();
    tr.col(2) = axes.col(2).normalized();
    tr_inv = tr.inverse();

    // Create array of local points in local basis
    MatrixXf coord(3,surf_spot.size());
    for(int j=0; j<surf_spot.size(); ++j){
        Vector3f p = sys->Box(frame).get_closest_image(surf_spot.XYZ(j),marker_coord);
        coord.col(j) = tr_inv*(p-marker_coord);
    }

    // Fit a quad surface
    Matrix<float,6,1> res;
    Vector3f sm; // Local smoothed coords
    fit_quad_surface(coord,res,sm,prop.fit_rms);

    // Get smoothed surface point in lab coords
    prop.smoothed = tr*sm + marker_coord;

    /* Now compute the curvatures

        First fundamental form:  I = E du^2 + 2F du dv + G dv^2
        E= r_u dot r_u, F= r_u dot r_v, G= r_v dot r_v

        For us parametric variables (u,v) are just (x,y) in local space.
        Derivatives:
            r_u = {1, 0, 2Ax+Cy+D}
            r_v = {0, 1, 2By+Cx+E}

        In central point x=0, y=0 so:
            r_u={1,0,D}
            r_v={0,1,E}

        Thus: E_ =1+D^2; F_ = D*E; G_ = 1+E^2;

        Second fundamental form: II = L du2 + 2M du dv + N dv2
        L = r_uu dot n, M = r_uv dot n, N = r_vv dot n

        Normal is just  n = {0, 0, 1}
        Derivatives:
            r_uu = {0, 0, 2A}
            r_uv = {0 ,0, C}
            r_vv = {0, 0, 2B}

        Thus: L_ = 2A; M_ = C; N_ = 2B;
        */

    float E_ = 1.0+res[3]*res[3];
    float F_ = res[3]*res[4];
    float G_ = 1.0+res[4]*res[4];

    float L_ = 2.0*res[0];
    float M_ = res[2];
    float N_ = 2.0*res[1];

    //Curvatures:
    prop.gaussian_curvature = (L_*N_-M_*M_)/(E_*G_-F_*F_);
    prop.mean_curvature = 0.5*(E_*N_-2.0*F_*M_+G_*L_)/(E_*G_-F_*F_);

    /* Compute principal directions and principal curvatures
         F1 = E F
              F G

         F2 = L M
              M N

         A = F1^-1*F2

         Eigenvalues of A are principal curvatures k1,k2, eigenvectors are principal directions
         */
    Matrix2f F1,F2;
    F1(0,0)=E_; F1(0,1)=F_;
    F1(1,0)=F_; F1(1,1)=G_;

    F2(0,0)=L_; F2(0,1)=M_;
    F2(1,0)=M_; F2(1,1)=N_;

    Eigen::SelfAdjointEigenSolver<Matrix2f> solver(F1.inverse()*F2);

    // We need to go back to lab coordinates for directions
    Matrix<float,3,2> lab_dirs;
    lab_dirs.block<2,2>(0,0) = solver.eigenvectors();
    lab_dirs(2,0) = lab_dirs(2,1) = 0.0; // local Z components
    lab_dirs.col(0).normalize();
    lab_dirs.col(1).normalize();
    lab_dirs = (tr*lab_dirs).eval(); // Go back to lab system

    prop.principal_curvatures = solver.eigenvalues();

    // If curvatures are negative the order of eigenvalues becomes inverted to what is needed
    // Thus order by absolute value
    int k1,k2;
    if(abs(prop.principal_curvatures(0))<abs(prop.principal_curvatures(1))){
        k2 = 0; k1 = 1;
    } else {
        k2 = 1; k1 = 0;
    }

    // Now k1 is along the largest curvature, k2 along the smallest

    prop.principal_directions.col(0) = lab_dirs.col(k1);
    prop.principal_directions.col(1) = lab_dirs.col(k2);

    return prop;
}

void Lipid_assembly::create(System &system, string lipids_sel, std::string heads_sel, string tails_sel, float d){
    System* sys = &system;

    Selection lipids(*sys,lipids_sel);

    // Select all heads
    Selection all_heads(*sys,heads_sel);
    // Sort into individual selections by resindex
    vector<Selection> heads;
    all_heads.split_by_residue(heads);

    // Select all tails
    Selection all_tails(*sys,tails_sel);
    // Sort into individual selections by resindex
    vector<Selection> tails;
    all_tails.split_by_residue(tails);

    // Get periodic centers of masses of all heads
    MatrixXf head_centers(3,heads.size());
    for(int i=0;i<heads.size();++i) head_centers.col(i) = heads[i].center(true,true);

    // We duplicate the frame and set first head atoms to centers of masses to get a
    // single point per lipid for grid searching
    sys->frame_dup(0);
    int fr = sys->num_frames()-1;
    Selection markers(*sys);
    for(int i=0;i<heads.size();++i) markers.append( heads[i].Index(0) );
    markers.set_frame(fr);
    for(int i=0;i<heads.size();++i){
        markers.XYZ(i,fr) = head_centers.col(i);
        heads[i].set_frame(fr);
        tails[i].set_frame(fr);
    }

    // Array of smoothed points
    MatrixXf smoothed(3,markers.size());

    ofstream vis_dirs("vis_local_dirs.tcl");

    markers.write("before_smooth.pdb");

    // Do several rounds of smoothing
    int Nsm = 10;
    for(int sm_iter=0; sm_iter<Nsm; ++sm_iter){

        cout << "Smoothing round #" << sm_iter << endl;

        // We need to re-wrap since centers of masses could be off box
        markers.wrap();

        // Find connectivity of head centers
        vector<Vector2i> bon;
        Grid_searcher(d,markers,bon,false,true);

        // Convert the list of bonds to convenient form
        vector<vector<int> > conn(heads.size());
        for(int i=0;i<bon.size();++i){
            conn[bon[i](0)].push_back(bon[i](1));
            conn[bon[i](1)].push_back(bon[i](0));
        }

        for(int i=0;i<markers.size();++i){
            // Selection for all neighbours including itself
            Selection local(*sys);
            local.append(markers.Index(i)); // append atom itself first
            for(int j=0; j<conn[i].size(); ++j){
                local.append(markers.Index(conn[i][j]));
            }
            local.set_frame(fr);

            Local_properties prop = get_local_curvature(local,tails[i]);
            smoothed.col(i) = prop.smoothed;

            // On last pass whow stats
            if(sm_iter==Nsm-1){
                cout << "gauss: " << prop.gaussian_curvature << endl;
                cout << "mean: " << prop.mean_curvature << endl;
                cout << "Principal curvatures: " << prop.principal_curvatures.transpose() << endl;
                cout << "Principal directions:\n" << prop.principal_directions.transpose() << endl;
                cout << "RMS: " << prop.fit_rms << endl;

                vis_dirs << tcl_arrow(markers.XYZ(i),markers.XYZ(i)+prop.principal_directions.col(0),0.5,"red");
                vis_dirs << tcl_arrow(markers.XYZ(i),markers.XYZ(i)+prop.principal_directions.col(1),0.5,"blue");
                vis_dirs << tcl_arrow(markers.XYZ(i),markers.XYZ(i)+prop.surface_normal,0.5,"green");

                markers.Occupancy(i) = prop.mean_curvature*100;
            }
        }

        // Reset positions to smoothed values
        markers.set_xyz(smoothed);

        markers.write("smoothed_"+boost::lexical_cast<string>(sm_iter)+".pdb");
    }

    vis_dirs.close();

}

