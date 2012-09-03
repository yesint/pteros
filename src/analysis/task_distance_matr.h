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

#ifndef TASK_DISTANCE_MATR_H
#define TASK_DISTANCE_MATR_H

#include "pteros/analysis/selection_analysis.h"

namespace pteros {

struct Task_distance_matr: public Task_base {
protected:
    Eigen::MatrixXf x1,x2,x3,x4; //accumulators for x and x**2
    int N;
    // For pairs within given distance
    float dist;
    Eigen::MatrixXf num;
    // Maximal moment
    int max_moment;
public:
    Task_distance_matr(Trajectory_processor* engine, Options_tree* opt);

    string task_name(){ return "distance_matr"; }
    int selections_required(){ return 1; }

    virtual void pre_process();
    virtual bool process_frame(const Frame_info& info);
    virtual void post_process(const Frame_info& info);
    //virtual void window_started_slot(const Trajectory_processor_info& info);
    //virtual void window_finished_slot(const Trajectory_processor_info& info);
    static void print_help();
};

}

#endif
