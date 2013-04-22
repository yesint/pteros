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

#include "pteros/analysis/trajectory_processor.h"
#include <string>
#include "pteros/analysis/message_channel.h"
#include "pteros/analysis/consumer.h"


using namespace std;
using namespace pteros;
using namespace Eigen;


Message_channel<std::string> ch;

void body1(){
    cout << "I'm thread " << boost::this_thread::get_id() << endl;
    for(int i=0;i<20;++i){
        cout << "Sending: Hi from thread! Iter: "+boost::to_string(i) << endl;
        ch.send("Hi from thread! Iter: "+boost::to_string(i));
    }
    ch.send_stop();
}

class Cons1: public Consumer {
public:
    Cons1(Trajectory_processor* pr): Consumer(pr){}
protected:
    virtual void pre_process(){
        cout << "(1) PRE " << endl;
    }
    virtual void post_process(const Frame_info& info){
        cout << "(1) POST " << info.absolute_frame << endl;
    }
    virtual bool process_frame(const Frame_info& info){
        cout << "(1) In frame " << info.absolute_frame << endl;
        return true;
    }
};

class Cons2: public Consumer {
public:
    Cons2(Trajectory_processor* pr): Consumer(pr){}
protected:
    virtual void pre_process(){
        cout << "(2) PRE " << endl;
    }
    virtual void post_process(const Frame_info& info){
        cout << "(2) POST " << info.absolute_frame << endl;
    }
    virtual bool process_frame(const Frame_info& info){
        cout << "(2) In frame " << info.absolute_frame << endl;
        boost::this_thread::sleep_for(boost::chrono::seconds(1));
        return true;
    }
};


int main(int argc, char** argv)
{
    /*
    try{
        Options_tree opt;
        opt.from_command_line(argc,argv);

        boost::thread worker(&body1);
        boost::this_thread::sleep_for(boost::chrono::seconds(2));
        System t("/home/semen/work/Projects/pteros/pteros_git/src/test/data/2lao.gro");
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


        //System t("/home/semen/work/Projects/pteros/pteros_git/src/test/data/2lao.gro");

    } catch(Pteros_error e){ e.print(); }
    */

    try{
        Options_tree opt;
        opt.from_command_line(argc,argv);
        Trajectory_processor proc(opt);
        Cons1 c1(&proc);
        //Cons2 c2(&proc);
        proc.run2();
    } catch(Pteros_error e){ e.print(); }

}

