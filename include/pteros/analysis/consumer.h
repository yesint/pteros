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


#ifndef CONSUMER_H
#define CONSUMER_H

#include <boost/signals2.hpp>
#include "pteros/core/system.h"
#include "pteros/simulation/simulation.h"
//#include "pteros/analysis/trajectory_processor.h"

namespace pteros {

/// Information about current frame, which is passed to Consumer for analysis
/// along with the frame itself
struct Frame_info {
    /// Counts *all* frames in trajectory group starting from 0
    int absolute_frame;
    /// Current time stamp
    int absolute_time;
    /// Counts only valid frames (the frames, which fall into specified time and frame range)
    /// and are sent to processing. Starts from 0.
    int valid_frame;
    /// Time of the first valid frame
    int first_time;
    /// Time of the last processed valid frame
    float last_time;
    /// First valid frame (this is absolute value!)
    int first_frame;
    /// Last processed valid frame (this is absolute value!)
    int last_frame;
    /// @name Window information
    /// {@
    /// Number of the current window
    int win_num;
    /// Size of window in time
    float win_size_time;
    /// Size of window in frames
    int win_size_frames;
    /// Valid frame when window started
    int win_start_frame;
    /// Time when window started
    float win_start_time;
    /// Valid frame when window finished (current if not finished yet)
    int win_last_frame;
    /// time when window finished (current if not finished yet)
    int win_last_time;
    /// @}
};

/// Auxiliary container for sending data to consumers
struct Data_container {
    /// Frame itself
    Frame frame;
    /// Frame information
    Frame_info frame_info;
    /// Stop flag. If it is true, than no frame is passed,
    /// but the end of processing is asked instead
    bool stop;
    Data_container(){
        stop = false;
    }
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
class Consumer {
public:
    Consumer(Trajectory_processor* pr);

    /// Main processing method. Called by Trajectory_processor
    virtual void run();
    System* get_system(){return &system;}
    void set_simulation(boost::shared_ptr<Simulation>& p){
        simulation = p;
    }

    //Simulation* get_simulation(){return &simulation;}
    void set_id(int i){id = i;}

    void set_proxy(bool pr){
        is_proxy = pr;
    }

    Frame* get_frame_ptr(){
        return proxy_frame_ptr;
    }

protected:
    /// This method is called before any processing starts and could be used
    /// to prepare anything for processing.
    /// This is a "housholding" method, which should not be overloaded by end user.
    virtual void setup();
    /// Called before process_frame for each frame
    /// This is a "housholding" method, which should not be overloaded by end user.
    virtual void before_each_frame();
    /// Called immediately before first frame is passed
    virtual void pre_process();
    /// Called immediately after last frame is processed
    virtual void post_process(const Frame_info& info);
    /// Called each time new frame arrives. This frame is stored in system.traj[0]
    virtual bool process_frame(const Frame_info& info);
    /// Called when new procesing window starts
    virtual void window_started(const Frame_info& info);
    /// Called when current procesing window ends
    virtual void window_finished(const Frame_info& info);

    /// Pointer to trajectory processor
    Trajectory_processor* proc;
    /// local system (stored in consumer itself)
    System system;
    /// Pointer to simulation (shared by all consumers)
    boost::shared_ptr<Simulation> simulation;
    /// Index of consumer
    int id;
    /// Window counter
    int win_num;

    void process_frame_info(Frame_info& info);
    float saved_time; //Save last processed timestamp
    int saved_abs_frame;
    int saved_valid_frame;

    bool is_proxy;
    Frame* proxy_frame_ptr;
};

}

#endif
