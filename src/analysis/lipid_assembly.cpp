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

using namespace std;
using namespace pteros;
using namespace Eigen;

Lipid_assembly::Lipid_assembly(Selection &sel, std::string head_marker_atom, float d){
    create(sel,head_marker_atom,d);
}

void Lipid_assembly::create(Selection &sel, std::string head_marker_atom, float d, float bilayer_cutoff){
    source_ptr = &sel;
    System* sys = sel.get_system();
    //Selection all(*sys,"all");
    // Create selection from marker atoms
    Selection markers(*sys,"("+sel.get_text()+") and "+head_marker_atom);    
    int N = markers.size();

    // Find connctivity of marker atoms
    vector<Vector2i> bon;
    Grid_searcher(d,markers,bon,false,true);

    // Convert the list of bonds to convenient form
    vector<vector<int> > conn(N);
    for(int i=0;i<bon.size();++i){
        conn[bon[i](0)].push_back(bon[i](1));
        conn[bon[i](1)].push_back(bon[i](0));
    }

    // Now cycle over connectivity map    
    Vector3f moments;
    Matrix3f axes;
    int atom;
    Selection local(*sys);          

    surface_normals.resize(markers.size());
    surface_mean_angle.resize(markers.size());
    surface_curvature.resize(markers.size());

    ofstream vis_dirs("vis_local_dirs.tcl");

    for(int i=0;i<markers.size();++i){
        cout << "(>>) " << i << endl;
        // Selection for all neighbours including itself
        local.clear();
        for(int j=0; j<conn[i].size(); ++j){
            local.append(markers.Index(conn[i][j]));
        }
        local.append(markers.Index(i)); // append atom itself

        // Compute the inertia axes        
        local.inertia(moments,axes,true);        
        // axes.col(2) is a local normal
        surface_normals[i] = axes.col(2);

        // Now fit a quadric surface to local points

        // First compute transformation matrix to local basis of inertia axes
        Matrix3f tr,tr_inv;
        tr.col(0) = axes.col(0).normalized();
        tr.col(1) = axes.col(1).normalized();
        tr.col(2) = axes.col(2).normalized();
        tr_inv = tr.inverse();

        // Create array of local points in local basis
        MatrixXf coord(3,local.size());
        for(int j=0; j<local.size(); ++j){
            Vector3f p = sys->get_closest_image(local.XYZ(j),markers.XYZ(i),0,false);
            coord.col(j) = tr_inv*(p-markers.XYZ(i));
        }

        //cout << coord.transpose() << endl;

        // We fit with polynomial fit = A*x^2 + B*y^2 + C*xy + D*x + E*y + F
        // Thus we need a linear system of size 6
        Matrix<float,6,6> m;
        Matrix<float,6,1> rhs,res; // Right hand side and result

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
                for(int j=0;j<local.size();++j){
                    m(r,c) += coef[r](j)*coef[c](j);
                }
            }
            // Now rhs
            for(int j=0;j<local.size();++j){
                rhs(r) += coef[r](j)*coord.col(j)(2);
            }
        }

        // Now solve
        res = m.colPivHouseholderQr().solve(rhs);

        // Compute RMS of fitting
        float rms = 0.0;
        for(int j=0;j<local.size();++j){
            float fit = 0.0;
            for(int r=0;r<6;++r) fit += coef[r](j)*res[r];
            rms += pow(coord.col(j)(2)-fit,2);
        }
        cout << "RMS: " << sqrt(rms/float(local.size())) << endl;

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
        float gaussian_curvature = (L_*N_-M_*M_)/(E_*G_-F_*F_);
        float mean_curvature = 0.5*(E_*N_-2.0*F_*M_+G_*L_)/(E_*G_-F_*F_);

        /* Compute principal directions and principal curvatures
         F1 = E F
              F G

         F2 = L M
              M N

         A = F1^-1*F2

         Eigenvalues of A are principal curvatures k1,k2, eigenvectors are principal directions
         */
        Matrix2f F1,F2;
        Vector2f princ_coef;
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


        princ_coef = solver.eigenvalues();

        cout << "gauss: " << gaussian_curvature << " mean: " << mean_curvature << endl;
        cout << "Principal curvatures: " << princ_coef.transpose() << endl;
        cout << "Principal directions:\n" << lab_dirs.transpose() << endl;

        vis_dirs << "draw color red" << endl;
        vis_dirs << "draw cylinder \"" << markers.XYZ(i).transpose()*10.0 << "\" ";
        vis_dirs << "\"" << (markers.XYZ(i)+lab_dirs.col(0)).transpose()*10.0 << "\" radius 0.5" << endl;

        vis_dirs << "draw color blue" << endl;
        vis_dirs << "draw cylinder \"" << markers.XYZ(i).transpose()*10.0 << "\" ";
        vis_dirs << "\"" << (markers.XYZ(i)+lab_dirs.col(1)).transpose()*10.0 << "\" radius 0.5" << endl;

        vis_dirs << "draw color green" << endl;
        vis_dirs << "draw cylinder \"" << markers.XYZ(i).transpose()*10.0 << "\" ";
        vis_dirs << "\"" << (markers.XYZ(i)+surface_normals[i]).transpose()*10.0 << "\" radius 0.5" << endl;


        markers.Occupancy(i) = abs(mean_curvature*100);
    }

    markers.write("colored.pdb");

    vis_dirs.close();

/*
    // Compute mean angle between the normals of local
    float a,ang;
    float c,curv;
    for(int i=0;i<markers.size();++i){
        //cout << "(**) " << i << endl;
        ang = 0;
        curv = 0;
        for(int j=0; j<conn[i].size(); ++j){
            // Angle between conn[i][j] and i
            a = acos(surface_normals[i].dot(surface_normals[conn[i][j]])/(surface_normals[i].norm() * surface_normals[conn[i][j]].norm()));
            if(a>M_PI_2) a = M_PI-a;
            ang += a;

            c = 2.0*sin(a) / sys->distance(markers.XYZ(i),markers.XYZ(conn[i][j]),0,true);
            curv += c;

        }
        if(conn[i].size()>0){
            ang /= conn[i].size();
            curv /= conn[i].size();
        }

        markers.Occupancy(i) = curv;

        surface_curvature[i] = curv;
        surface_mean_angle[i] = ang;
    }

    // All markers with occupancy < cutoff correspond to bilayer part
    // Obtain monolayers
    Selection bi(*sys,"occupancy < " + boost::lexical_cast<string>(bilayer_cutoff));
    markers.write("colored.pdb");
    */
}

