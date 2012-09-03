/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
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

#include "task_box.h"
#include <fstream>

using namespace std;
using namespace pteros;

void Task_box::pre_process(){
    data.clear();
    volume.clear();
}

bool Task_box::process_frame(const Frame_info &info){
    Eigen::Matrix3f m = system.Box(0);
    data.push_back(m.diagonal());
    volume.push_back(m.diagonal().prod());
    cout << "***" << info.absolute_time << " " << Selection(system,"all").XYZ(0).transpose() << endl;
    return true;
}

void Task_box::post_process(const Frame_info &info){
    // Output
    string fname = prefix+".dat";
    // Get time step in frames and time
    float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);
    //cout << "BOX: " << info.last_time << " " << info.first_time << " " << info.valid_frame << endl;

    ofstream f(fname.c_str());
    f << "# time box_a box_b box_c box_volume:" << endl;
    for(int i=0; i<data.size(); ++i){
        f << i*dt << " " << data[i].transpose() << " " << volume[i] << endl;
    }
    f.close();
}

//void Task_box::window_started_slot(const Trajectory_processor_info& info){}
//void Task_box::window_finished_slot(const Trajectory_processor_info& info){}

void Task_box::print_help(){
    cout << "Task box:\n"
            "-----------------\n"
            "Computes volume of the rectangular periodic box and a,b,c box vectors.\n"
         << endl;
}
