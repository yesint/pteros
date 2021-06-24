/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#include "pteros/extras/gnm.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "pteros/core/distance_search.h"
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <ctime>

using namespace std;
using namespace pteros;
using namespace Eigen;

GNM::GNM(const Selection &sel, float cutoff){
    int i,j,k;

    N = sel.size();
    eigenvalues.resize(N-1);
    eigenvectors.resize(N,N-1);

    LOG()->debug("Constructing GNM Kirkgoff matrix. Cut-off: {}, Size: {}",cutoff,N);

    MatrixXf kirk(N,N); //Kirgoff matrix
    float d, time1, time2;

    kirk.fill(0.0);

    vector<Eigen::Vector2i> bon;
    vector<float> dist;
    search_contacts(cutoff,sel,bon,dist,false,noPBC);

    for(int i=0; i<bon.size(); ++i){
        kirk(bon[i](0),bon[i](1)) = kirk(bon[i](1),bon[i](0)) = -1.0;
    }

    // Compute diagonal elements
    for(i=0;i<N;++i) kirk(i,i) = -kirk.col(i).sum();

    LOG()->debug("Computing eigenvectors...");

    time1 = clock();
    Eigen::SelfAdjointEigenSolver<MatrixXf> solver(kirk);
    eigenvalues = solver.eigenvalues();
    eigenvectors = solver.eigenvectors();
    // Vectors are already sorted, cool :)
    time2 = clock();

    LOG()->debug(" Done in {} s.",(time2-time1)/CLOCKS_PER_SEC);
}

Eigen::VectorXf GNM::get_eigenvector(int i) const
{
    return eigenvectors.col(i);
}

VectorXf GNM::get_B_factor() const
{
    return b;
}

void GNM::write_eigenvectors(string fname, int v1, int v2){
    if(v1<0 || v2>N-2 || v2<v1) throw PterosError("Can't write these eigenvectors!");
    ofstream f(fname.c_str());
    if(!f) throw PterosError("Can't open file "+fname+" for writing!");
    int i,j;
    for(i=0;i<N;++i){
        f << i << " ";
        for(j=v1;j<=v2;++j) f << eigenvectors(i,j) << " ";
        f << endl;
    }
    f.close();
}

void GNM::compute_c_matrix(){
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

    // Save B-factors
    b.resize(N);
    b = c.diagonal();

    // Normalize by self-correlations to get correct cross-correlations
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
    if(!f) throw PterosError("Can't open file "+fname+" for writing!");
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
    if(!f) throw PterosError("Can't open file "+fname+" for writing!");

    for(i=0;i<N;++i){
        for(j=0;j<N;++j)
            f << p(i,j) << " ";
        f << endl;
    }
    f.close();
}

MatrixXf GNM::get_subset_c_matrix(ArrayXi_const_ref subset) const
{
    int sz = subset.size();
    MatrixXf ret(sz,sz);

    for(int i=0;i<sz;++i){
        for(int j=0;j<sz;++j){
            ret(i,j) = c(subset[i],subset[j]);
        }
    }

    return ret;
}

MatrixXf GNM::get_c_matrix() const
{
    return c;
}
