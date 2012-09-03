#include "task_covar_matr.h"
#include <fstream>
#include <set>

using namespace std;
/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
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

using namespace pteros;
using namespace Eigen;


Task_covar_matr::Task_covar_matr(Trajectory_processor *engine, Options_tree *opt):
    Task_base(engine,opt)
{    
    do_align = options->get_value<bool>("align",true);
}

void Task_covar_matr::pre_process(){
    system.frame_dup(0);

    N = sel[0].size();

    matr.resize(3*N,3*N); matr.fill(0.0);
    mean.resize(3,N); mean.fill(0.0);    
}

bool Task_covar_matr::process_frame(const Frame_info &info){
    // Align if requested
    if(do_align) sel[0].fit(0,1);

    int i,j;
    // Do accumulation
    for(i=0;i<N;++i){
        mean.col(i) += sel[0].XYZ(i);
        for(j=i;j<N;++j){
            matr.block<3,3>(3*i,3*j) += sel[0].XYZ(i) * sel[0].XYZ(j).transpose();
        }
    }

    return true;
}

void Task_covar_matr::post_process(const Frame_info &info){
    string fname;
    ofstream f;

    float T = float(info.valid_frame+1);

    // Get average
    mean /= T;
    matr.triangularView<Eigen::Upper>() /= T;

    int i,j;
    // Now compute covariance
    for(i=0;i<N;++i){
        for(j=i;j<N;++j){
            matr.block<3,3>(3*i,3*j) -= mean.col(i)*mean.col(j).transpose();
        }
    }

    // Make simmetric
    matr.triangularView<Eigen::StrictlyLower>() = matr.triangularView<Eigen::StrictlyUpper>().transpose();

    bool write_whole = options->get_value<bool>("write_whole",false);
    if(write_whole){
        // Output whole matrix
        fname = prefix+"-whole.dat";
        f.open(fname.c_str());
        f << "# Covariance matrix (3N*3N)" << endl;
        f << matr << endl;
        f.close();
    }

    // Output matrix of diagonal elements only (cov(x)+cov(y)+cov(z))
    Eigen::MatrixXd m(N,N);
    for(i=0;i<N;++i){
        for(j=0;j<N;++j){
            m(i,j) = matr.block<3,3>(3*i,3*j).trace();
        }
    }

    fname = prefix+".dat";
    f.open(fname.c_str());
    f << "# Covariance matrix (N*N)" << endl;
    f << m << endl;
    f.close();
}

//void Task_box::window_started_slot(const Trajectory_processor_info& info){}
//void Task_box::window_finished_slot(const Trajectory_processor_info& info){}

void Task_covar_matr::print_help(){
    cout << "Task covar_matr:\n"
            "-----------------\n"
            "Computes covariance matrix for selection.\n"
         << endl;
}
