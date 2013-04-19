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

#include "pteros/pteros.h"
#include "pteros/core/grid_search.h"
#include "pteros/analysis/trajectory_processor.h"
#include "pteros/analysis/rmsf.h"
#include "pteros/analysis/bilayer.h"
#include "pteros/core/mol_file.h"
#include <string>

#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include "pteros/analysis/message_channel.h"


using namespace std;
using namespace pteros;
using namespace Eigen;

Message_channel<std::string> ch;

void body1(){
    cout << "I'm thread " << std::this_thread::get_id() << endl;
    for(int i=0;i<20;++i){
        cout << "Sending: Hi from thread! Iter: "+boost::to_string(i) << endl;
        ch.send("Hi from thread! Iter: "+boost::to_string(i));
    }
    ch.send_stop();
}


int main(int argc, char** argv)
{
    try{        
        std::thread worker(body1);
        std::this_thread::sleep_for(std::chrono::seconds(2));
        //System t("/home/semen/work/Projects/pteros/pteros_git/src/test/data/2lao.gro");
        //System t("/media/data/semen/trajectories/test/topol.tpr.pttop");
        //t.load("/media/data/semen/trajectories/test/traj.trr");
        string s;

        while(ch.recieve(s)){
            cout << "Got message: " << s << endl;
        }

        while(!ch.empty()){
            ch.recieve(s);
            cout << "Got message: " << s << endl;
        }

        worker.join();

    } catch(Pteros_error e){ e.print(); }

}

