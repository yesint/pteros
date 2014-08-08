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

#include "pteros/core/gnm.h"
#include "pteros/core/pteros_error.h"
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <ctime>

using namespace std;
using namespace pteros;
using namespace Eigen;

void GNM::compute(Selection& sel, float cutoff){
    int i,j,k;

    cout << "Running GNM..."<<endl;

    N = sel.size();
    eigenvalues.resize(N-1);
    eigenvectors.resize(N,N-1);

    cout << "Constructing GNM Kirkgoff matrix. Cut-off: " << cutoff << " Size: " << N << endl;
    MatrixXf kirk(N,N); //Kirgoff matrix
    float d, time1, time2;

    kirk.fill(0.0);
    // Compute off-diagonal elements
    for(i=0;i<N-1;++i)        
        for(j=i+1;j<N;++j){            
            d = (sel.XYZ(i)-sel.XYZ(j)).norm();
            if(d<=cutoff){
                kirk(i,j)= -1.0;
                kirk(j,i)= -1.0;
            }
        }

cout << "**" << endl;
    // Compute diagonal elements
    for(i=0;i<N;++i) kirk(i,i) = -kirk.col(i).sum();

    cout << "Computing eigenvectors...";
    time1 = clock();

    Eigen::SelfAdjointEigenSolver<MatrixXf> solver(kirk);
    eigenvalues = solver.eigenvalues();
    eigenvectors = solver.eigenvectors();
    // Vectors are already sorted, cool :)

    time2 = clock();
    cout << " Done in " << (time2-time1)/CLOCKS_PER_SEC << " s." << endl;
}

GNM::GNM(Selection& sel, float cutoff){
    compute(sel,cutoff);
}

void GNM::write_eigenvectors(string fname, int v1, int v2){
    if(v1<0 || v2>N-2 || v2<v1) throw Pteros_error("Can't write these eigenvectors!");
    ofstream f(fname.c_str());
    if(!f) throw Pteros_error("Can't open file "+fname+" for writing!");
    int i,j;
    for(i=0;i<N;++i){
        f << i << " ";
        for(j=v1;j<=v2;++j) f << eigenvectors(i,j) << " ";
        f << endl;
    }
    f.close();
}

void GNM::compute_c_matrix(bool normalize){
    int i,j,k;
    cout << "Calculating correlations..." << endl;
    c.resize(N,N);
    c.fill(0.0);
    for(k=1;k<N;++k) // For all eigenvectors except first, which is zero
        for(i=0;i<N;++i)
            for(j=i;j<N;++j)
                c(i,j) = c(i,j) + eigenvectors(i,k)*eigenvectors(j,k)/eigenvalues(k);

    // Set lower triangle of the matrix
    for(i=0;i<N-1;++i)
        for(j=i+1;j<N;++j)
            c(j,i) = c(i,j);

    if(normalize)
        for(i=0;i<N;++i)
            for(j=0;j<N;++j)
                c(i,j) = c(i,j)/sqrt(c(i,i)*c(j,j));
}

void GNM::compute_p_matrix(){
    if(!c.size()) compute_c_matrix();
    cout << "Calculating correlations of correlation patterns...";
    int i,j;
    float time1,time2;

    time1 = clock();

    p.resize(N,N);
    p.fill(0.0);

    // See if c-matrix is ready
    if(!c.cols()) compute_c_matrix();

    // Pre-compute m and s
    VectorXf m(N), s(N);

    for(i=0;i<N;++i){
        m(i) = c.col(i).sum()/float(N);
        s(i) = c.col(i).array().pow(2).sum()/float(N) - pow(m(i),2);
    }

    for(i=0;i<N;++i)
        for(j=i;j<N;++j){
            p(i,j) = (c.col(i).array() * c.col(j).array()).sum() /float(N) - m(i)*m(j);
            p(i,j) = p(i,j) / sqrt(s(i)*s(j));
            p(j,i) = p(i,j);
        }

    time2 = clock();
    cout << " Done in " << (time2-time1)/CLOCKS_PER_SEC << " s." << endl;
}

void GNM::write_c_matrix(string fname){
    // See if c-matrix is ready
    if(!c.cols()) compute_c_matrix();

    int i,j;
    ofstream f(fname.c_str());
    if(!f) throw Pteros_error("Can't open file "+fname+" for writing!");
    for(i=0;i<N;++i){
        for(j=0;j<N;++j)
            f << c(i,j) << " ";
        f << endl;
    }
    f.close();
}

void GNM::write_p_matrix(string fname){
    // See if p-matrix is ready
    if(!p.cols()) compute_p_matrix();

    int i,j;
    ofstream f(fname.c_str());
    if(!f) throw Pteros_error("Can't open file "+fname+" for writing!");

    for(i=0;i<N;++i){
        for(j=0;j<N;++j)
            f << p(i,j) << " ";
        f << endl;
    }
    f.close();
}
