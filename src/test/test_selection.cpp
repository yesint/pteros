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

#include "pteros/core/peg_parser.h"

#include <list>

using namespace std;
using namespace pteros;
using namespace Eigen;



//--------------------------------------------------------------------------------------

int main(int argc, char** argv)
{

    try{        
/*
        Parser_data INP;
        //INP.string_to_parse = "-( 2.3e-8+ (x - 3.014 ) *1.2)-(3+4- -4*5.5-(2+-3.1)) ";
        INP.string_to_parse = "name CA";
        INP.cur_buf_pos = 0;
        GREG g;        
        yyinit(&g);        
        g.data = &INP;
        while (yyparse(&g));
        if(g.begin != g.end) throw Pteros_error("Syntax error at pos "+to_string(g.begin));

        yydeinit(&g);
        INP.root->dump();
*/
        //System s("/home/semen/work/Projects/asymmetric_bilayer/for-diamonds/hexagonal/2x2.gro");
        System s1("/home/semen/work/Projects/asymmetric_bilayer/for-diamonds/hexagonal/128.pdb");
        Selection sel(s1,"(name CA CB 'N.*' and index 1 2 3 5-25) or within 2.3 of (x^2+y^2>2.5+43.2*(y+z)/3)");


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

