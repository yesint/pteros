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

using namespace std;
using namespace pteros;
using namespace Eigen;


int main(int argc, char** argv)
{

    try{
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
        */

        Options toplevel;
        vector<Options> tasks;

        parse_command_line(argc,argv,toplevel,"task",tasks);
        //parse_command_line(argc,argv,toplevel);
        toplevel.print();
        cout << "------" << endl;
        for(auto& t: tasks){
            t.print();
            cout << "------" << endl;
        }

        float v = toplevel["t"].as_float();
        cout << v << endl;


    } catch(const Pteros_error& e){ e.print(); }

}

