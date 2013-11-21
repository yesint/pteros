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
#include <stack>
#include <sstream>

using namespace std;
using namespace pteros;
using namespace Eigen;

Lipid_assembly::Lipid_assembly(System &system, std::string lipids_sel, string heads_sel, string tails_sel, float d, int Nsmooth){
    create(system,lipids_sel,heads_sel,tails_sel,d, Nsmooth);
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

float angle_between_vectors(const Ref<const Vector3f>& v1, const Ref<const Vector3f>& v2){
    return acos(v1.dot(v2)/(v1.norm() * v2.norm()));
}

void Lipid_assembly::get_local_curvature(Selection& surf_spot, Selection* tail_spot, Local_curvature& prop, bool force_flip){

    System* sys = surf_spot.get_system();
    int frame = surf_spot.get_frame();

    // Central marker atom (by convension the first atom in surf_spot)
    Vector3f marker_coord = surf_spot.XYZ(0);

    // Compute the inertia axes
    Vector3f moments;
    Matrix3f axes;
    surf_spot.inertia(moments,axes,true);

    // If tails selection is provided orient normals according to it
    if(tail_spot){
        // Compute the center of masses of tails ends
        Vector3f end = sys->Box(frame).get_closest_image(tail_spot->center(true,true),marker_coord);
        Vector3f lip_vec = marker_coord-end;

        // Force reversal is asked for
        if(force_flip) lip_vec = -lip_vec;

        // Now an angle between the lipid vector and normal
        float a = angle_between_vectors(axes.col(2).normalized(),lip_vec);
        //float a = acos(prop.surface_normal.dot(lip_vec)/(prop.surface_normal.norm() * lip_vec.norm()));
        if(a>M_PI_2){ // Need to flip normal
            axes.col(2) = -axes.col(2).eval();
        }
    }

    // axes.col(2) is a local normal
    prop.surface_normal = axes.col(2).normalized();

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
    prop.point = tr*sm + marker_coord;

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
}

void Lipid_assembly::create_markers(vector<Selection>& lipids,Selection& markers){
    // Get periodic centers of masses        
    MatrixXf lip_centers(3,lipids.size());
    for(int i=0;i<lipids.size();++i) lip_centers.col(i) = lipids[i].center(true,true);

    // We duplicate the frame and set first head atoms to centers of masses to get a
    // single point per lipid for grid searching
    lipids[0].get_system()->frame_dup(markers.get_frame());
    ++extra_frames;   

    int fr = lipids[0].get_system()->num_frames()-1;
    markers.modify(*lipids[0].get_system());

    for(int i=0;i<lipids.size();++i) markers.append( lipids[i].Index(0), false );

    markers.set_frame(fr);
    for(int i=0;i<lipids.size();++i){
        markers.XYZ(i,fr) = lip_centers.col(i);
    }

}

void Lipid_assembly::write_vmd_arrows(Selection& markers, const vector<Local_curvature>& props, string fname){
    ofstream vis_dirs(fname);
    for(int i=0;i<markers.size();++i){
        vis_dirs << tcl_arrow(markers.XYZ(i),markers.XYZ(i)+props[i].principal_directions.col(0),0.5,"red");
        vis_dirs << tcl_arrow(markers.XYZ(i),markers.XYZ(i)+props[i].principal_directions.col(1),0.5,"blue");
        vis_dirs << tcl_arrow(markers.XYZ(i),markers.XYZ(i)+props[i].surface_normal,0.5,"green");
    }
    vis_dirs.close();
}

void Lipid_assembly::write_mean_curvature_as_pdb(Selection& markers, const vector<Local_curvature>& props, string fname){
    for(int i=0;i<markers.size();++i){
        markers.Occupancy(i) = props[i].mean_curvature*100;
    }
    markers.write(fname);
}

Local_curvature Lipid_assembly::weighted_curvature_in_point(Vector3f &point,
                                                            float cutoff,
                                                            Selection *spot)
{
    Local_curvature curv;
    curv.point = point;

    System* sys = head_markers.get_system();
    int fr_h = head_markers.get_frame();
    int fr_t = tail_markers.get_frame();
    // Find lipids around this point
    set<int> m;
    vector<int> bon;

    float cur_dist;

    cur_dist = (cutoff<0) ? dist : cutoff;

    do {
        Grid_searcher g;
        g.assign_to_grid(cur_dist,head_markers,false,true);
        g.search_within(point,bon);

        if(bon.size()==0){
            //cout << "(WARNING!) no head markers within " << cur_dist << "! Increasing to " << cur_dist*1.5 << endl;
            cur_dist *= 1.5;
        }
    } while(bon.size()==0);

    // Make spot selection if asked
    if(spot){
        spot->modify(*sys);
        for(int i=0; i<bon.size(); ++i) spot->append(head_markers.Index(bon[i]),false);
    }

    // Compare orientations of the normals
    // First point serves as a reference
    VectorXi normal_dir(bon.size());
    normal_dir(0) = 1; // Reference is up
    float a;
    int n_up = 1, n_down = 0;
    for(int n=1; n<bon.size(); ++n){
        a = angle_between_vectors(head_props[bon[n]].surface_normal, head_props[bon[0]].surface_normal);
        if(a>M_PI_2){
            normal_dir(n) = -1; // Opposite to reference
            ++n_down;
        } else {
            normal_dir(n) = 1; // Same as reference
            ++n_up;
        }
    }
    // If down is more common flip everything
    if(n_down>n_up) normal_dir = -normal_dir;

    float sum_dist = 0.0, max_dist = -1.0;
    VectorXf dist_h(bon.size()),dist_t(bon.size());

    for(int n=0; n<bon.size(); ++n){
        //dist_h(n) = sys->Box(fr_h).distance(point,head_markers.XYZ(bon[n]));
        //dist_t(n) = sys->Box(fr_t).distance(point,tail_markers.XYZ(bon[n]));
        dist_h(n) = sys->Box(fr_h).distance(point,head_props[bon[n]].point);
        dist_t(n) = sys->Box(fr_t).distance(point,tail_props[bon[n]].point);
        if(dist_h(n)>max_dist) max_dist = dist_h(n);
        if(dist_t(n)>max_dist) max_dist = dist_t(n);
    }

    // Distance from the up surface is returned in roughness field!

    for(int n=0; n<bon.size(); ++n){
        // Head of lipid
        sum_dist += max_dist-dist_h(n);
        curv.gaussian_curvature += head_props[bon[n]].gaussian_curvature * (max_dist-dist_h(n)) * normal_dir(n);
        curv.mean_curvature += head_props[bon[n]].mean_curvature * (max_dist-dist_h(n)) * normal_dir(n);
        curv.surface_normal += head_props[bon[n]].surface_normal * (max_dist-dist_h(n)) * normal_dir(n);

        // Tail of lipid
        sum_dist += max_dist-dist_t(n);
        curv.gaussian_curvature += tail_props[bon[n]].gaussian_curvature * (max_dist-dist_t(n)) * normal_dir(n);
        curv.mean_curvature += tail_props[bon[n]].mean_curvature * (max_dist-dist_t(n)) * normal_dir(n);
        curv.surface_normal += tail_props[bon[n]].surface_normal * (max_dist-dist_t(n)) * normal_dir(n);
    }

    curv.gaussian_curvature /= sum_dist;
    curv.mean_curvature /= sum_dist;
    curv.surface_normal = (curv.surface_normal/sum_dist).normalized().eval();

    // Find minimal distance to up point
    float min_dist = 1e20;
    int min_ind = -1;
    for(int n=0; n<bon.size(); ++n){
        if(normal_dir(n)>0){ // Only count up surface for distance!
            if(dist_h(n)<min_dist){
                min_dist = dist_h(n);
                min_ind = n;
            }
        }
    }
    // vector from point to corresponding head
    Vector3f v = head_props[bon[min_ind]].point -
            sys->Box(fr_h).get_closest_image(point,head_props[bon[min_ind]].point);
    // Project v onto weighted normal --> distance
    curv.roughness = v.dot(curv.surface_normal);

    //cout << tcl_arrow(point,point+curv.surface_normal,0.5,"yellow") << endl;

    return curv;
}

void Lipid_assembly::compute_surface(Selection& markers, vector<Selection>& tails,
                                     float d, int Nsm,
                                     vector<Local_curvature>& props, bool force_flip){
    //markers.write("before_smooth.pdb");
    System* sys = markers.get_system();
    int fr = markers.get_frame();

    // Array of smoothed points
    MatrixXf smoothed(3, markers.size());

    // Do several rounds of smoothing
    for(int sm_iter=0; sm_iter<Nsm; ++sm_iter){

        cout << "\tIteration " << sm_iter << " of " << Nsm-1 << endl;

        // We need to re-wrap since centers of masses could be off box
        markers.wrap();

        // Find connectivity of markers
        vector<Vector2i> bon;
        Grid_searcher(d,markers,bon,false,true);

        // Convert the list of bonds to convenient form
        vector<vector<int> > conn(markers.size());
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

            get_local_curvature(local,&tails[i],props[i],force_flip);
            smoothed.col(i) = props[i].point;

            // On first pass save roughness
            if(sm_iter==0){
                props[i].roughness = props[i].fit_rms;
            }

        }

        // Reset positions to smoothed values
        markers.set_xyz(smoothed);

        //markers.write("smoothed_"+boost::lexical_cast<string>(sm_iter)+".pdb");
    }
}


void Lipid_assembly::create(System &system, string lipids_sel, std::string heads_sel, string tails_sel, float d, int Nsmooth){
    System* sys = &system;

    dist = d;
    Nsm = Nsmooth;

    lipids.modify(*sys,lipids_sel);

    // Select all heads
    Selection all_heads(*sys,heads_sel);       
    all_heads.split_by_residue(heads);

    // Select all tails
    Selection all_tails(*sys,tails_sel);  
    all_tails.split_by_residue(tails);

    head_markers.modify(*sys);
    tail_markers.modify(*sys);

    head_props.resize(heads.size());
    tail_props.resize(heads.size());

    extra_frames = 0;
}

void Lipid_assembly::compute(int frame){

    clock_t t1 = clock();
    // If this is not the first invocation than there are two
    // extra frames created on last invocation
    // Delete them:
    if(extra_frames>0){
        cout << "Deleting " << extra_frames << " old aux frames" << endl;
        head_markers.get_system()->frame_delete(head_markers.get_system()->num_frames()-extra_frames);
        extra_frames = 0;
    }

    lipids.set_frame(frame);
    for(int i=0;i<heads.size();++i) heads[i].set_frame(frame);
    for(int i=0;i<heads.size();++i) tails[i].set_frame(frame);
    head_markers.set_frame(frame);
    tail_markers.set_frame(frame);

    cout << "Smoothing heads..." << endl;
    create_markers(heads,head_markers);
    compute_surface(head_markers,tails,dist,Nsm,head_props);

    cout << "Smoothing tails..." << endl;
    create_markers(tails,tail_markers);
    // Tails surface should have the same normal as corresponding heads
    // so we need to force reverse of the normals
    compute_surface(tail_markers,heads,dist,Nsm,tail_props,true);    

    clock_t t2 = clock();
    cout << "Computation time: " << (t2-t1)/(float)CLOCKS_PER_SEC << endl;
}

void Lipid_assembly::write_output(){
    write_vmd_arrows(head_markers,head_props,"vis_local_dirs_heads.tcl");
    write_mean_curvature_as_pdb(head_markers,head_props,"heads_smoothed.pdb");

    write_vmd_arrows(tail_markers,tail_props,"vis_local_dirs_tails.tcl");
    write_mean_curvature_as_pdb(tail_markers,tail_props,"tails_smoothed.pdb");
}
