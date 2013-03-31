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
#include "pteros/simulation/simulation.h"
#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

namespace pteros {

/** The base class for trajectory processing
*   It provides facilities for loading large trajectories by chunks
*   and to analyze each frame by user-defined function.
*   The range of processing could be given
*   by frame number or by physical time.
*/

typedef boost::shared_ptr<Frame> Frame_ptr;

class Trajectory_processor {
    friend class Consumer_base;
    public:

        /// Default constructor
        Trajectory_processor(){ init(); }
        /// Constructor with system
        Trajectory_processor(Options_tree& opt){
            init();
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
        static void print_help();

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
        /// Simulation object, which is global for all consumers
        boost::shared_ptr<Simulation> simulation;

        // Sets default options
        void init();           

        int log_interval;

        /// Async reading stuff
        boost::thread reading_thread; //reading thread itself
        boost::thread_group consumer_threads; // Consumer threads
        int max_buffer_size; // Maximal buffer size
        // Frame buffer. Contains shared pointers to data containers!
        std::map<int,boost::shared_ptr<Data_container> > buffer;
        boost::mutex buffer_mutex; // Buffer mutex
        boost::condition buffer_cond;
        bool stop_requested;        
        std::map<int,std::vector<bool> > frame_access_count; // Who accessed this frame?

        void read_single_trajectory(std::string& fname);
        int abs_frame; // Absolte frame index
        int valid_frame; // Valid frame index
        // Saved first frame and time
        int saved_first_frame;
        float saved_first_time;
        bool check_time_range(int fr, float t);

        std::vector<Consumer_base*> consumers;
        std::vector<bool> alive_consumers;

        void consumer_finished(int id);
        boost::shared_ptr<Data_container> frame_provider(int fr, int id);
        void add_consumer(Consumer_base* p);
        void start_threads(std::vector<std::string>& fnames);

        /// Reader parameters
        int first_frame, last_frame;
        float first_time, last_time;
        int skip;
        int window_size_frames;
        float window_size_time;
        float custom_start_time;
        float custom_dt;
        void fill_window_info(Frame_info& info);
};

}
#endif

