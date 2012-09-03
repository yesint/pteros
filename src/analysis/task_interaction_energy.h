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

#ifndef TASK_INTERACTION_ENERGY_H
#define TASK_INTERACTION_ENERGY_H

#include "pteros/analysis/selection_analysis.h"

namespace pteros {

struct Task_interaction_energy: public Task_base {
    std::vector<Energy_components> data;
    Energy_components mean;
    float cutoff;
    // Grid searcher
    Grid_searcher searcher;
    bool is_periodic;
    std::vector<Eigen::Vector2i> bon;

    Task_interaction_energy(Trajectory_processor* engine, Options_tree* opt);

    string task_name(){ return "interaction_energy"; }
    int selections_required(){ return 2; }

    virtual void pre_process();
    virtual bool process_frame(const Frame_info& info);
    virtual void post_process(const Frame_info& info);

protected:
    bool is_self_energy;
};

}

#endif
