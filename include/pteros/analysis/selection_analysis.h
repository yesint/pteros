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


#ifndef SELECTION_ANALYSIS_H
#define SELECTION_ANALYSIS_H

#include "pteros/analysis/trajectory_processor.h"
#include "pteros/core/selection.h"
#include "pteros/simulation/simulation.h"
#include <Eigen/Core>
#include <boost/ptr_container/ptr_vector.hpp>

namespace pteros {

/// Base class for all analysis tasks. Subclass of Consumer.
struct Task_base: public Consumer {
    /** Constructor.
      @arg proc Trajectory_processor, which will be bound to this task
      @arg opt Options, which are passed to this task and could be interpreted inside
      @arg set_texts Selection strings of all selections, which are used by this task
      @arg sel_name Symbolic names of all selections
      */
    Task_base(Trajectory_processor* proc, Options_tree* opt);
    void create(std::vector<std::string>& sel_texts,
                std::vector<std::string>& sel_names);

    /// Returns name of the task (should be defined in derived classes!)
    virtual std::string task_name() = 0;

    /// Returns number of required selections (-1 means all supplied)
    virtual int selections_required() = 0;

    bool task_ready(){ return is_ready; }

protected:
    Options_tree* options;

    /// Overrided setup code
    virtual void setup();

    /// Symbolic names of selections
    std::vector<std::string> sel_name;

    /// Selection texts
    std::vector<std::string> sel_text;

    /// Selections
    std::vector<Selection> sel;

    /// Unique prefixes for output files
    std::string prefix;    

    /// Override to update selections on each frame
    virtual void before_each_frame();

    bool is_ready;
};


/// Driver, which runs multiple tasks and executes them in parallel
class Selection_analysis
{
public:
    Selection_analysis();
    Selection_analysis(Trajectory_processor& proc, Options_tree& opt){
        create(proc,opt);
    }

    void create(Trajectory_processor& proc, Options_tree& opt);
    static void print_help();

    /// Factory for creating tasks
    boost::shared_ptr<Task_base> task_factory(std::string task_name, Options_tree* opt);
private:
    Trajectory_processor* engine;
    Options_tree* options;

    /// List of tasks
    vector<boost::shared_ptr<Task_base> > tasks;
};




}

#endif // SELECTION_ANALYSIS_H


