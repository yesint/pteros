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

#include <string>
#include "pteros/analysis/options.h"
#include <Eigen/Core>
#include "pteros/core/pteros_error.h"
#include "pteros/core/selection.h"
#include "pteros/core/distance_search.h"

#include "pteros/analysis/trajectory_reader.h"
#include "pteros/analysis/task_base.h"

#include <chrono>
#include <fstream>
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;

class SelTest_mask {
public:
    SelTest_mask(const Selection& sel){
        mask.resize(sel.get_system()->num_atoms());
        mask.fill(0);
        n = sel.size();
        for(int i=0;i<sel.size();++i){
            mask(sel.Index(i)) = 1;
        }
        b = sel.Index(0);
        e = sel.Index(sel.size()-1);
        coord = &sel.get_system()->Frame_data(sel.get_frame()).coord;
    }

    ~SelTest_mask(){}

    Vector3f center(){
        Vector3f res(Vector3f::Zero());
        for(int i=b; i<=e; ++i)
            if(mask(i)==1)
                res += (*coord)[i];
        return res/float(n);
    }

private:
    VectorXi mask;
    int b,e,n;
    vector<Vector3f>* coord;
};


//-------------------------------------------------------

TASK_PARALLEL(Test_task)
    virtual void pre_process(){
        cout << "Test_task pre_process" << endl;
    }
    virtual void process_frame(const Frame_info& info){

        cout << "Test_task process_frame " << std::this_thread::get_id() << " "<< info.valid_frame << endl;
        //std::this_thread::sleep_for(std::chrono::seconds(1));
        for(int i=0; i<10; ++i)
            system().rotate(1,0.02);
    }

    virtual void post_process(const Frame_info& info){
        cout << "Test_task post_process of instance " << info.last_frame << endl;
    }
};

void accum(const Frame_info& info, const std::vector<Task_ptr>& tasks){
    cout << "Running collector for " << tasks.size() << " tasks" << endl;
}


TASK_SERIAL(Test_serial)
    virtual void pre_process(){
        cout << "Test_serial pre_process" << endl;
    }
    virtual void process_frame(const Frame_info& info){

        cout << "Test_serial process_frame " << std::this_thread::get_id() << " "<< info.valid_frame << endl;
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    virtual void post_process(const Frame_info& info){
        cout << "Test_serial post_process" << info.last_frame << endl;
    }
};



int main(int argc, char** argv)
{

    try{        

        Options opt;
        parse_command_line(argc,argv,opt);
        Trajectory_reader reader(opt);


        reader.add_task( new Test_task() );
        reader.register_collector( &accum );

        //reader.add_task( new Test_serial() );
        //reader.add_task( new Test_serial() );

        reader.run();

        //System s("/home/semen/work/Projects/Besancon-2014/Graphene/gr.pdb");
        //vector<Vector2i> bon;
        //Selection sel(s,"all");

    /*
        vector<float> dist;
        search_contacts(0.143,sel,bon,true,true,&dist);

        for(int i=0;i<bon.size();++i){
            cout << bon[i].transpose() << " " << dist[i] << endl;
        }
*/
        //cout << (boost::get<Parse_tree_ptr>(p->children.front())) << endl;

        //std::shared_ptr<Parser> p(new Parser);
        //p->parse();

//        System s("/home/semen/work/Projects/asymmetric_bilayer/for-diamonds/hexagonal/2x2.gro");
  //      Selection sel(s,"(y+4)<(x+2) or (1-z)>x");
/*
        System s("/home/semen/work/Projects/pteros/paper2/supplement/benchmark/large.gro");


        int Ntry = 1000;


        Vector3f r1(Vector3f::Zero()), r2(Vector3f::Zero());

        ofstream f("size-dep.dat");

        for(int N=2;N<s.num_atoms();N+=50000){
            // Make sparse vector
            vector<int> ind;
            ind.reserve(N);
            for(int i=0;i<N;++i){
                ind.push_back(i*(s.num_atoms()/(N-1)));
            }

            // Normal selection
            r1.fill(0);
            Selection sel(s,ind);

            auto t_start = std::chrono::high_resolution_clock::now();
            for(int n=0;n<Ntry;++n)
                r1 += sel.center();
            auto t_end = std::chrono::high_resolution_clock::now();
            auto T1 = 1e6*std::chrono::duration<double>(t_end-t_start).count()/float(Ntry);

            // Mask access
            r2.fill(0);
            SelTest_mask sel1(sel);

            t_start = std::chrono::high_resolution_clock::now();
            for(int n=0;n<Ntry;++n)
                r2 += sel1.center();
            t_end = std::chrono::high_resolution_clock::now();
            auto T2 = 1e6*std::chrono::duration<double>(t_end-t_start).count()/float(Ntry);

            cout << N << " " << T1 << " " << T2 << " | " << r1.transpose() << " :: " << r2.transpose() << endl;
            f << N << " " << T1 << " " << T2 << endl;
        }

        f.close();
*/

/*
        Selection w;
        Selection sel(s,"all");
        int N = 10;
        vector<Vector2i> bon;
        float d = 0.6;

        auto t_start = std::chrono::high_resolution_clock::now();
        for(int i=0;i<N;++i)
            w.modify(s,"within 0.6 nopbc of name PO4");
            //Grid_searcher(d,sel,bon);
        auto t_end = std::chrono::high_resolution_clock::now();

        cout << " elapsed: "
             << 1e6*std::chrono::duration<double>(t_end-t_start).count()/float(N) << endl;
*/

    } catch(const Pteros_error& e){ e.print(); }

}

