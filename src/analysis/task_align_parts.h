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

#ifndef TASK_ALIGN_PARTS_H
#define TASK_ALIGN_PARTS_H

#include "pteros/analysis/selection_analysis.h"
#include "pteros/core/mol_file.h"

namespace pteros {

struct Task_align_parts: public Task_base {
private:
    std::string out_fname;
    boost::shared_ptr<Mol_file> trj;
public:
    Task_align_parts(Trajectory_processor* engine, Options_tree* opt);

    string task_name(){ return "align_parts"; }
    int selections_required(){ return -1; }

protected:
    virtual void pre_process();
    virtual bool process_frame(const Frame_info& info);
    virtual void post_process(const Frame_info& info);

    Selection all;
    //virtual void window_started_slot(const Trajectory_processor_info& info);
    //virtual void window_finished_slot(const Trajectory_processor_info& info);
};

}

#endif
