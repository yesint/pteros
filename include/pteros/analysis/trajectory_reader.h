/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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


#ifndef TRAJECTORY_PROCESSOR_H
#define TRAJECTORY_PROCESSOR_H

#include <string>
#include <memory>
#include <functional>
#include "pteros/core/selection.h"
#include "pteros/analysis/options.h"
#include "pteros/analysis/consumer_base.h"

namespace pteros {

typedef std::shared_ptr<Frame> Frame_ptr;
typedef Message_channel<std::shared_ptr<Data_container> > Data_channel;
typedef std::shared_ptr<Data_channel> Data_channel_ptr;



class Task_base {
    friend class Trajectory_reader;
public:
    Task_base(){}
    virtual ~Task_base(){}
    virtual Task_base* clone() const = 0;
protected:
    virtual bool is_parallel() = 0;
    virtual void pre_process() = 0;
    virtual void process_frame(const Frame_info& info) = 0;
    virtual void post_process(const Frame_info& info) = 0;

private:

    void put_frame(const Frame& frame){
        system.Frame_data(0) = frame;
    }

    void put_system(const System& sys){
        system = sys;
    }

    System system;
};


typedef std::shared_ptr<Task_base> Task_ptr;

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
        virtual ~Trajectory_reader(){};

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

        int log_interval;

        bool is_frame_valid(int fr, float t);
        bool is_end_of_interval(int fr, float t);

        /// Reader parameters
        int first_frame, last_frame;
        float first_time, last_time;
        int skip;

        float custom_start_time;
        float custom_dt;

        void reader_thread_body(const Data_channel_ptr &channel);
        std::vector<std::string> traj_files;

        void process_value_with_suffix(const std::string& s, int* intval, float* floatval);

        std::vector<Task_ptr> tasks;

        bool is_parallel;

        Frame_info dispatch_frames_to_task(const Task_ptr &task,
                                            const Data_channel_ptr &channel,
                                            const System& sys);

        std::function<void(const Frame_info&,const std::vector<Task_ptr>&)> collector;

};

}
#endif

