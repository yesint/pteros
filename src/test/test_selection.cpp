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
#include "pteros/core/grid_search.h"
#include <chrono>

using namespace std;
using namespace pteros;
using namespace Eigen;


int main(int argc, char** argv)
{

    try{

        //System s("/media/semen/data/semen/trajectories/asymmetric_hexagonal/with_c60/last.gro");
        System s("/media/semen/data/semen/trajectories/2lao/after_em.gro");
        vector<Vector2i> bon;

        /*
        auto t_start = std::chrono::high_resolution_clock::now();
        Grid_searcher(2.0,sel1,sel2,bon,true,true);
        auto t_end = std::chrono::high_resolution_clock::now();

        cout << bon.size() << " elapsed: "
             << std::chrono::duration<double>(t_end-t_start).count() << endl;
        */
        //-----------

        Selection sel(s,"all");
        int N=100000;

        Vector3f v(1,2,3);

        VectorXi mask(s.num_atoms());
        mask.fill(0);
        for(int i=0;i<sel.size();++i) mask(sel.Index(i))=1;



        auto t_start = std::chrono::high_resolution_clock::now();
        for(int i=0;i<N;++i){
            sel.translate(v);
        }
        auto t_end = std::chrono::high_resolution_clock::now();
        cout << " elapsed: "
             << std::chrono::duration<double>(t_end-t_start).count()/double(N) << endl;

        /*
        auto t_start = std::chrono::high_resolution_clock::now();
        Selection w;
        for(int i=0;i<100;++i)
            w.modify(s,"within 4.0 noself nopbc of name CA");
        auto t_end = std::chrono::high_resolution_clock::now();

        cout << w.size() << " elapsed: "
             << std::chrono::duration<double>(t_end-t_start).count()/100.0 << endl;
        */


        /*
        Selection sel(s,"not name CA");
        Selection sel2(s,"name CA");
        Grid_searcher(1.0,sel,sel2,bon,true,false);
        cout << bon.size() << endl;
        */

        /*
        tng_trajectory_t trj;
        string fname("/home/semen/work/Projects/Besancon-2014/cisplatin/fit/amber/traj.tng");
        tng_util_trajectory_open(fname.c_str(),'r',&trj);

        int stat = TNG_SUCCESS;
        char datatype;
        void* values = 0;
        int64_t frame_num;
        double frame_time;

        while(stat==TNG_SUCCESS){

            stat = tng_util_particle_data_next_frame_read(trj, TNG_TRAJ_POSITIONS, &values,
                                                          &datatype, &frame_num, &frame_time);
            cout << stat << " " << frame_num<< endl;
        }
        */

        //System s("/media/semen/data/semen/trajectories/asymmetric_hexagonal/with_c60/last.gro");
        //System s("/home/semen/work/Projects/Besancon-2014/cisplatin/fit/amber/after_md.pdb");
        //s.load("/home/semen/work/Projects/Besancon-2014/cisplatin/fit/amber/traj.tng");
        //System s("/home/semen/work/Projects/Besancon-2014/cisplatin/fit/amber/traj.tng");

        //Selection sel;
        //sel.modify(s,"x>1");
        //cout << sel.size() << endl;
        //sel.write("test.tng");

/*
Grammar g("ups");
g.root();

while(!g.actions.empty()){
    g.actions.front()();
    g.actions.remove(g.actions.front());
}
*/
        //Selection sel(s,"name CA");
/*
std::unique_ptr<T> p,p1;
p.reset(new T);
p1.reset(new T);
string s("name CA");
p->doit(s);
s = "name BF";
p1->doit(s);
*/
        /*
        string str("--trajectory["
                   "initial_structure.pdb "
                   "traj.xtc r1/traj.xtc r1/traj.part0002.xtc r1/traj.part0003.xtc "
                   "--first_frame 0"
                   "--last_frame 100 "
                   "--log_interval 2 "
                   "] "
                   "--task rms[--selection \"name CA\" --unwrap 0.2] "
                   "--task rms[--selection \"protein\" --unwrap 0.4] "
                   "--dump_input dump");
        cout << str << endl;

        Options_tree opt;
        opt.from_command_line(str);
        for(auto o: opt.get_options("task")){
            cout << o->get_value<string>("selection") << endl;
            cout << o->get_value<double>("unwrap") << endl;
        }

        System s;


        Options toplevel;
        vector<Options> tasks;

        parse_command_line(argc,argv,toplevel,"task",tasks);
        //parse_command_line(argc,argv,toplevel);
        toplevel.debug();
        cout << "------" << endl;
        for(auto& t: tasks){
            t.debug();
            cout << "------" << endl;
        }

        vector<float> v = toplevel("tramvay","3.14159 42 -4.5").as_floats();
        v = toplevel("tramvay").as_floats();
        for(auto a: v) cout << a << endl;

*/
    } catch(const Pteros_error& e){ e.print(); }

}

