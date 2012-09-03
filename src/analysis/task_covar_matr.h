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

#ifndef TASK_COVAR_MATR_H
#define TASK_COVAR_MATR_H

#include "pteros/analysis/selection_analysis.h"

namespace pteros {

struct Task_covar_matr: public Task_base {
protected:
    Eigen::MatrixXf matr,mean;
    int N;
    bool do_align;

public:
    Task_covar_matr(Trajectory_processor* engine, Options_tree* opt);

    string task_name(){ return "covar_matr"; }
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
