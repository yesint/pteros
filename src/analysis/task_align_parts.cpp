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

#include "task_align_parts.h"
#include <fstream>

using namespace std;
using namespace pteros;


Task_align_parts::Task_align_parts(Trajectory_processor *engine, Options_tree *opt):
    Task_base(engine,opt)
{    
    out_fname = options->get_value<string>("output_file",prefix+"out_traj.xtc");    
}

void Task_align_parts::pre_process(){
    cout << "Will align " << sel.size() << " selections " << endl;
    system.frame_dup(0);
    trj = io_factory(out_fname,'w');
    if(!trj->get_content_type().trajectory)
        throw Pteros_error("File "+out_fname+" is not a trajectory!");

    all.modify(system,"all");
}

bool Task_align_parts::process_frame(const Frame_info &info){
    Eigen::Affine3f t;
    for(int i=0;i < sel.size(); ++i){        
        sel[i].fit(0,1);
    }
    Mol_file_content c;
    c.trajectory = true;
    trj->write(all,c);
    return true;
}

void Task_align_parts::post_process(const Frame_info &info){
    trj.reset(); // Will release Mol_file writer
}

//void Task_center::window_started_slot(const Trajectory_processor_info& info){}
//void Task_center::window_finished_slot(const Trajectory_processor_info& info){}

