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

#include "task_center.h"
#include <fstream>

using namespace std;
using namespace pteros;


Task_center::Task_center(Trajectory_processor *engine, Options_tree *opt):
    Task_base(engine,opt)
{    
    weighted = options->get_value<bool>("weighted",false);
}

void Task_center::pre_process(){
    cout << "in preprocess\n";
    mean.fill(0.0);
    data.clear();
}

bool Task_center::process_frame(const Frame_info &info){
    sel[0].set_frame(0);
    Eigen::Vector3f v = sel[0].center(weighted);
    data.push_back(v);
    mean += v;
    return true;
}

void Task_center::post_process(const Frame_info &info){
    mean /= (float)info.valid_frame;
    // Output
    string fname = prefix+".dat";
    // Get time step in frames and time
    float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);
    //cout << "CENTER: " << info.last_time << " " << info.first_time << " " << info.valid_frame << endl;

    ofstream f(fname.c_str());
    f << "# Center of selection '" << sel_name[0] << "' [" << sel_text[0] << "]" << endl;
    f << "# Mean: " << mean.transpose() << endl;
    f << "# time X Y Z:" << endl;
    for(int i=0; i<data.size(); ++i){
        f << i*dt << " " << data[i].transpose() << endl;
    }
    f.close();
}

//void Task_center::window_started_slot(const Trajectory_processor_info& info){}
//void Task_center::window_finished_slot(const Trajectory_processor_info& info){}

