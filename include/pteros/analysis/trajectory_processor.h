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

#include <iostream>
#include <string>
#include <queue>
#include <set>
#include "pteros/core/selection.h"
#include <fstream>

#include "pteros/analysis/options_parser.h"
#include "pteros/analysis/consumer_base.h"
#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

#include "pteros/analysis/message_channel.h"

namespace pteros {

/** The base class for trajectory processing
*   It provides facilities for loading large trajectories by frames
*   and to analyze each frame by user-defined function.
*   The range of processing could be given
*   by frame number or by physical time.
*/

typedef std::shared_ptr<Frame> Frame_ptr;
typedef Message_channel<std::shared_ptr<Data_container> > Data_channel;

class Trajectory_processor {
    friend class Consumer_base;
    public:

        /// Default constructor
        Trajectory_processor(){ }
        /// Constructor with system
        Trajectory_processor(Options_tree& opt){            
            options = &opt;
        }
        /// Destructor
        virtual ~Trajectory_processor();

        /// Pass options
        void set_options(Options_tree& opt){
            options = &opt;
        }

        /// Do computation
        virtual void run();

        /// Print summary of allowed options    
        std::string help();

    protected:
        /// Tree of options
        /** Tree of options is supposed to have the following structure:
            root:
                log_interval:
                    int
                trajectory_group:
                    range: (optional)
                        type:
                            "frame_range" | "time_range"
                        begin:
                            int | double
                        end:
                            int | double
                    trajectory_file_1
                    trajectory_file_2
                    trajectory_file_3...
                trajectory_group:
                    ...
        */
        Options_tree* options;

        int log_interval;

        bool is_frame_valid(int fr, float t);
        bool is_end_of_interval(int fr, float t);

        std::vector<Consumer_base*> consumers;

        void add_consumer(Consumer_base* p);        

        /// Reader parameters
        int first_frame, last_frame;
        float first_time, last_time;
        int skip;
        int window_size_frames;
        float window_size_time;
        float custom_start_time;
        float custom_dt;

        Data_channel channel;

        void reader_thread_body();
        std::vector<std::string> traj_files;
};

}
#endif

