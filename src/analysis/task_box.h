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

#ifndef TASK_BOX_H
#define TASK_BOX_H

#include "pteros/analysis/selection_analysis.h"

namespace pteros {

/// Task for computing dimensions and volume of the periodic box.
/// Empty list of selections should be passed, because box is not selection-dependent.
struct Task_box: public Task_base {
    std::vector<Eigen::Vector3f> data;
    std::vector<float> volume;

    Task_box(Trajectory_processor* engine, Options_tree* opt): Task_base(engine,opt) {}

    string task_name(){ return "box"; }
    int selections_required(){ return 0; }

    virtual void pre_process();
    virtual bool process_frame(const Frame_info& info);
    virtual void post_process(const Frame_info& info);
    //virtual void window_started_slot(const Trajectory_processor_info& info);
    //virtual void window_finished_slot(const Trajectory_processor_info& info);
    static void print_help();
};

}

#endif
