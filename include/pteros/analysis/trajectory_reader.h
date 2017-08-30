/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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


#ifndef TRAJECTORY_READER_H
#define TRAJECTORY_READER_H

#include <string>
#include <functional>
#include "pteros/analysis/options.h"
#include "pteros/analysis/task_base.h"


namespace pteros {

typedef std::shared_ptr<Task_base> Task_ptr;
typedef Message_channel<std::shared_ptr<pteros::Data_container> > Data_channel;
typedef std::shared_ptr<Data_channel> Data_channel_ptr;

/** The base class for trajectory processing
*   It provides facilities for loading large trajectories by frames
*   and to analyze each frame by user-defined function.
*   The range of processing could be given
*   by frame number or by physical time.
*/
class Trajectory_reader {
public:

        /// Default constructor
        Trajectory_reader();
        /// Constructor with options
        Trajectory_reader(const Options& opt);
        /// Destructor
        virtual ~Trajectory_reader(){}

        /// Pass options
        void set_options(const Options& opt){
            options = opt;
        }

        /// Read trajectory
        virtual void run();

        /// Print summary of allowed options
        std::string help();

        /// Adds new task
        void add_task(Task_base* task){
            tasks.push_back( Task_ptr(task) );
        }

        /// Register collecting function for parallel tasks
        void register_collector( std::function<void(const Frame_info&,const std::vector<Task_ptr>&)> func){
            collector = func;
        }

private:

        // Options
        Options options;

        //void reader_thread_body(const Data_channel_ptr &channel);

        std::vector<std::string> traj_files;        

        std::vector<Task_ptr> tasks;

        bool is_parallel;

        std::function<void(const Frame_info&,const std::vector<Task_ptr>&)> collector;

};

}
#endif

