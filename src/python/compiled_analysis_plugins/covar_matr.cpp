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


#include "pteros/python/compiled_plugin.h"
#include <fstream>

using namespace std;
using namespace pteros;
using namespace Eigen;

class covar_matr: public pteros::Compiled_plugin_base {
public:
    covar_matr(pteros::Trajectory_processor* pr, const pteros::Options& opt): Compiled_plugin_base(pr,opt) {}

protected:

    void pre_process(){
        do_align = options("align","true").as_bool();        
        sel.modify(system,options("selection").as_string() );

        system.frame_dup(0);

        N = sel.size();

        matr.resize(3*N,3*N); matr.fill(0.0);
        mean.resize(3,N); mean.fill(0.0);
    }

    void process_frame(const Frame_info &info){
        // Align if requested
        if(do_align) sel.fit(0,1);

        int i,j;
        // Do accumulation
        for(i=0;i<N;++i){
            mean.col(i) += sel.XYZ(i);
            for(j=i;j<N;++j){
                matr.block<3,3>(3*i,3*j) += sel.XYZ(i) * sel.XYZ(j).transpose();
            }
        }
    }

    void post_process(const Frame_info &info){
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

        bool write_whole = options("write_whole","false").as_bool();
        if(write_whole){
            // Output whole matrix
            fname = label+"-whole.dat";
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

        fname = label+".dat";
        f.open(fname.c_str());
        f << "# Covariance matrix (N*N)" << endl;
        f << m << endl;
        f.close();
    }

private:
    bool do_align;
    int N;
    MatrixXf matr;
    MatrixXf mean;
    Selection sel;
};

CREATE_COMPILED_PLUGIN(covar_matr)




