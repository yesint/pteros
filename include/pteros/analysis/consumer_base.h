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

#ifndef CONSUMER_BASE_H
#define CONSUMER_BASE_H

#include "pteros/core/system.h"
#include "pteros/analysis/frame_info.h"
#include "pteros/analysis/message_channel.h"

#include <fstream>

#include <iostream>

namespace pteros {

/// Auxiliary container for sending data to consumers
struct Data_container {
    /// Frame itself
    Frame frame;
    /// Frame information
    Frame_info frame_info;
};

// Forward declaration
class Trajectory_processor;

/** Base class for asynchronous data analysis.
  Consumer recieves frames one by one and process them by means
  of virtual functions, which should be overloaded in subclasses.
  Creating consumer only makes sence in conjunction with Trajectory_processor,
  which provides frame for processing.
  When Consumer is created it is registered in Trajectory_processor automatically.
  */
class Consumer_base {
    friend class Trajectory_processor;
public:
    Consumer_base(Trajectory_processor* pr);

    System* get_system(){return &system;}
    void set_id(int i){id = i;}    

protected:
    /// Called immediately before first frame is passed
    virtual void pre_process();
    /// Called immediately after last frame is processed
    virtual void post_process(const Frame_info& info);
    /// Called each time new frame arrives. This frame is stored in system.traj[0]
    virtual void process_frame(const Frame_info& info);
    /// Called when new procesing window starts
    virtual void window_started(const Frame_info& info);
    /// Called when current procesing window ends
    virtual void window_finished(const Frame_info& info);

    /// local system (stored in consumer itself)
    System system;
    /// Pointer to trajectory processor
    Trajectory_processor* proc;

    /// Removes jumps on each frame for given set of atoms,
    /// which is stored in no_jump_ind selection
    /// If no_jumps_ind is empty (default) does nothing
    /// User have to populate it in pre_process() callback using
    /// add_no_jump_atoms() methods
    void remove_jumps(const Frame_info &info);

    /// Handler functions which call user callbacks
    /// Could be overriden to take additional actions
    /// These handlers are called by Trajectory_processor
    virtual void pre_process_handler(){
        in_pre_process = true;
        in_process_frame = in_post_process = false;
        pre_process();
    }

    virtual void post_process_handler(const Frame_info& info){
        in_post_process = true;
        in_pre_process = in_process_frame = false;
        post_process(info);
    }

    virtual void process_frame_handler(const Frame_info& info){
        in_process_frame = true;
        in_pre_process = in_post_process = false;        
        // Call user callback
        process_frame(info);
    }

    virtual void window_started_handler(const Frame_info& info){
        window_started(info);
    }

    virtual void window_finished_handler(const Frame_info& info){
        window_finished(info);
    }

    // Flags, which indicate where user-called method is executed
    bool in_pre_process, in_process_frame, in_post_process;

private:

    /// Index of consumer
    int id;
    /// Window counter
    int win_num;
    float win_start_time;
    int win_start_frame;

    void process_window_info(Frame_info& info);

    // Should be defined in derived classes
    virtual void process_frame_data(Frame& data) = 0;

    float saved_time; //Save last processed timestamp
    int saved_abs_frame;
    int saved_valid_frame;    

    void run_in_thread(std::shared_ptr<Message_channel<std::shared_ptr<Data_container> > >& chan);
    void consume_frame(std::shared_ptr<Data_container>& data);    
};

}

#endif
