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

#include "task_rmsd.h"
#include <fstream>

using namespace std;
using namespace pteros;


Task_rmsd::Task_rmsd(Trajectory_processor *engine, Options_tree *opt):
    Task_base(engine,opt)
{    
}

void Task_rmsd::pre_process(){
    mean = 0.0;
    data.clear();
}

bool Task_rmsd::process_frame(const Frame_info &info){
    // Fitting breaks the system, but we have local copy, nobody cares. Cool :)
    // Set reference frame
    if(info.valid_frame==0){
        system.frame_dup(0);
    }

    Eigen::Affine3f trans = sel[0].fit_transform(0,1);
    sel[0].apply_transform(trans);
    float v = sel[0].rmsd(0,1);
    data.push_back(v);
    mean += v;

    return true;
}

void Task_rmsd::post_process(const Frame_info &info){
    mean /= (float)info.valid_frame;
    // Output
    string fname = prefix+".dat";
    // Get time step in frames and time
    float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);

    ofstream f(fname.c_str());
    f << "# RMSD of selection '" << sel_name[0] << "' [" << sel_text[0] << "]" << endl;
    f << "# Mean: " << mean << endl;
    f << "# time X Y Z:" << endl;
    for(int i=0; i<data.size(); ++i){
        f << i*dt << " " << data[i] << endl;
    }
    f.close();
}

//void Task_rmsd::window_started_slot(const Trajectory_processor_info& info){}
//void Task_rmsd::window_finished_slot(const Trajectory_processor_info& info){}

