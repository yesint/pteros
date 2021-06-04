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
#include "pteros/core/distance_search.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

class distance_matr: public pteros::Compiled_plugin_base {
public:
    distance_matr(pteros::Trajectory_processor* pr, const pteros::Options& opt): Compiled_plugin_base(pr,opt) {}

protected:

    void pre_process(){
        // See if computation of pairs within given distance is requested (0 means not requested)
        dist = options("dist","0.0").as_float();
        max_moment = options("max_moment","4.0").as_int();
        if(max_moment<1 || max_moment>4) throw Pteros_error("Moments from 1 to 4 are supported");        
        sel.modify(system,options("selection").as_string() );

        system.frame_dup(0);

        N = sel.size();
        if(max_moment>=1){
            x1.resize(N,N);
            x1.fill(0.0);
        }
        if(max_moment>=2){
            x2.resize(N,N);
            x2.fill(0.0);
        }
        if(max_moment>=3){
            x3.resize(N,N);
            x3.fill(0.0);
        }
        if(max_moment==4){
            x4.resize(N,N);
            x4.fill(0.0);
        }

        if(dist>0){
            num.resize(N,N);
            num.fill(0);
        }
    }

    void process_frame(const Frame_info &info){

        // Do accumulation
        int i,j;
        float r;
        if(dist>0){
            vector<Eigen::Vector2i> bon;
            search_contacts(dist,sel,bon);
            for(i=0;i<bon.size();++i){
                if(abs(bon[i](0)-bon[i](1))>4){
                    r = (sel.XYZ(bon[i](0))-sel.XYZ(bon[i](1))).norm();
                    if(max_moment>=1) x1(bon[i](0),bon[i](1)) += r;
                    if(max_moment>=2) x2(bon[i](0),bon[i](1)) += r*r;
                    if(max_moment>=3) x3(bon[i](0),bon[i](1)) += r*r*r;
                    if(max_moment==4) x4(bon[i](0),bon[i](1)) += r*r*r*r;
                    ++num(bon[i](0),bon[i](1));
                }
            }
        } else {
            Eigen::Vector3f p;
            for(i=0;i<N-1;++i){
                p = sel.XYZ(i);
                for(j=i+1;j<N;++j){
                    r = (sel.XYZ(j)-p).norm();
                    if(max_moment>=1) x1(i,j) += r;
                    if(max_moment>=2) x2(i,j) += r*r;
                    if(max_moment>=3) x3(i,j) += r*r*r;
                    if(max_moment==4) x4(i,j) += r*r*r*r;
                }
            }
        }
    }

    void post_process(const Frame_info &info){
        string fname;
        ofstream f;

        float T = float(info.valid_frame+1);

        if(dist>0){
            // Compute mean and dispersion
            if(max_moment>=1)
                x1 = (num.array() > 0).select(x1.array()/num.array(), 0.0);
            if(max_moment>=2)
                x2 = (num.array() > 0).select(x2.array()/num.array() - x1.array()*x1.array(), 0.0);
            if(max_moment>=3)
                x3 = (num.array() > 0).select(
                        2.0*x1.array()*x1.array()*x1.array() - 3.0*x1.array()*x2.array() + x3.array()/T
                  ,0.0);
            if(max_moment==4)
                x4 = (num.array() > 0).select(
                    -3.0*x1.array().pow(4) + 6.0*x1.array().pow(2)*x2.array()
                    -4.0*x1.array()*x3.array() + x4.array()/T
                  ,0.0);

        } else {
            // Finish computing mean. Place it to x
            if(max_moment>=1)
                x1 = x1.array()/T;
            // Finish computing variance. Place it to x2
            if(max_moment>=2)
                x2 = x2.array()/T - x1.array()*x1.array();
            // Compute mu3
            if(max_moment>=3)
                x3 = 2.0*x1.array()*x1.array()*x1.array() - 3.0*x1.array()*x2.array() + x3.array()/T;
            // Compute mu4
            if(max_moment>=4)
                x4 = -3.0*x1.array().pow(4) + 6.0*x1.array().pow(2)*x2.array()
                    -4.0*x1.array()*x3.array() + x4.array()/T;
        }

        if(max_moment>=3){
            // Compute skewness and place to x3
            x3 = (x2.array() != 0).select(x3.array()/x2.array().pow(3/2), 0.0);
        }
        if(max_moment>=4){
            // Compute kurtosis and place to x4
            x4 = (x2.array() != 0).select(x4.array()/x2.array().pow(2)-3.0, 0.0);
        }

        if(max_moment>=1) x1.triangularView<Eigen::StrictlyLower>() = x1.triangularView<Eigen::StrictlyUpper>().transpose();
        if(max_moment>=2) x2.triangularView<Eigen::StrictlyLower>() = x2.triangularView<Eigen::StrictlyUpper>().transpose();
        if(max_moment>=3) x3.triangularView<Eigen::StrictlyLower>() = x3.triangularView<Eigen::StrictlyUpper>().transpose();
        if(max_moment==4) x4.triangularView<Eigen::StrictlyLower>() = x4.triangularView<Eigen::StrictlyUpper>().transpose();

        // Output
        fname = label+"-mean.dat";
        f.open(fname.c_str());
        f << "# Mean distance matrix" << endl;
        f << x1 << endl;
        f.close();

        if(max_moment>=2){
            fname = label+"-disp.dat";
            f.open(fname.c_str());
            f << "# Distance dispersion matrix" << endl;
            f << x2 << endl;
            f.close();
        }

        if(max_moment>=3){
            fname = label+"-skew.dat";
            f.open(fname.c_str());
            f << "# Distance skewness matrix" << endl;
            f << x3 << endl;
            f.close();
        }

        /*
        // Skewness profile
        VectorXf v = x3.colwise().sum();
        v = 100.0*(v.array()-v.minCoeff())/(v.maxCoeff()-v.minCoeff());
        for(int i=0; i<N; ++i) sel[0].Beta(i) = v(i);
        sel[0].write("colored-skewness.pdb");
        */
        if(max_moment==4){
            fname = label+"-kurtosis.dat";
            f.open(fname.c_str());
            f << "# Distance kurtosis matrix" << endl;
            f << x4 << endl;
            f.close();
        }
    }

private:
    float dist;
    int max_moment, N;
    Selection sel;
    MatrixXf x1,x2,x3,x4,num;
};

CREATE_COMPILED_PLUGIN(distance_matr)




