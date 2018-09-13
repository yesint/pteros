/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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


#include <string>
#include "pteros/analysis/options.h"
#include <Eigen/Core>
#include "pteros/core/pteros_error.h"
#include "pteros/core/selection.h"
#include "pteros/core/distance_search.h"


#include "pteros/analysis/trajectory_reader.h"
#include "pteros/analysis/task_plugin.h"

#include <chrono>
#include <fstream>
#include <thread>

#include "pteros/core/logging.h"
#include "spdlog/fmt/ostr.h"

#include "pteros/extras/membrane.h"
#include "pteros/extras/substructure_search.h"

#include <ctime>

using namespace std;
using namespace pteros;
using namespace Eigen;

//-----------------------------

int main(int argc, char** argv)
{

    /*
    string t1("name C21 C22 C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216");
    string t2("name C21 C22 C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216");

    vector<Lipid_descr> species = {
        {"YOPE","resname YOPE", "name P N", "name C314 C315 C316 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"DYPE","resname DYPE", "name P N", "name C314 C315 C316 C214 C215 C216", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"DOPE","resname DOPE", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"DOPC","resname DOPC", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"POPE","resname POPE", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"DYPC","resname DYPC", "name P N", "name C314 C315 C316 C214 C215 C216", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"YOPC","resname YOPC", "name P N", "name C314 C315 C316 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"POPC","resname POPC", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
    };

    System sys("/home/semen/work/current/Projects/Masato/symmetric/after_em.gro");

    Membrane membr(&sys,species);
    membr.compute_properties(2.0, true, Vector3f(0,0,0), Vector3i(0,0,1));
    membr.write_vmd_arrows("normals.tcl");
    membr.write_smoothed("smoothed.pdb");

    // lip properties

    for(int i=0;i<membr.lipids.size();++i){
        //cout << Map<VectorXf>(membr.lipids[i].order[0].data(),16).transpose() << endl;
        cout << membr.lipids[i].area << endl;
    }

    return 1;
    */

    //System s("/home/semen/work/current/Projects/plasma_vesicle/after_em.gro");
/*
    // System with 10^6 atoms
    System s("/home/semen/work/stored/Projects/pteros/paper2/supplement/benchmark/large.gro");
    vector<int> ind;
    for(int i=0;i<s.num_atoms();i+=1) ind.push_back(i);
    Selection sel(s,ind);

    std::clock_t start;
    double duration;

    //-----------------------

    {
        start = std::clock();

        vector<Vector3f> data(sel.size());
        for(int i=0;i<sel.size();++i) data[i] = sel.xyz(i);


        for(int i=1;i<sel.size();++i){
            for(int j=1;j<1000;++j){
                data[i]=data[i-1]+data[j];
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;
    }

    //-----------------------

    {
        start = std::clock();

        MatrixXf data(3,sel.size());
        for(int i=0;i<sel.size();++i) data.col(i) = sel.xyz(i);

        for(int i=1;i<sel.size();++i){
            for(int j=1;j<1000;++j){
                data.col(i)=data.col(i-1)+data.col(j);
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;
    }

    //-----------------------
    {

        start = std::clock();

        for(int i=1;i<sel.size();++i){
            for(int j=1;j<1000;++j){
                sel.xyz(i)=sel.xyz(i-1)+sel.xyz(j);
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;
    }

    //-----------------------


    {

        start = std::clock();

        for(int i=1;i<sel.size();++i){
            for(int j=1;j<1000;++j){
                sel[i].xyz()=sel[i-1].xyz()+sel[j].xyz();
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;
    }
    //-----------------------

    return 1;
*/
    try{


        System src("/home/semen/work/current/Projects/Ache/b.pdb");
        System sample("/home/semen/work/current/Projects/Ache/b_sample.pdb");

        auto res = find_substructures(src(),sample(),true);
        cout << res.size() << endl;

        //src_sel.write("/home/semen/work/current/Projects/Ache/1.mol2");
        //Topmatch t(src_sel);


    } catch(const Pteros_error& e){
        LOG()->error(e.what());
    }

}


