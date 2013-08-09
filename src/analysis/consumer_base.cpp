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


#include "pteros/analysis/consumer_base.h"
#include "pteros/analysis/trajectory_processor.h"
#include <boost/bind.hpp>

using namespace std;
using namespace pteros;

Consumer_base::Consumer_base(Trajectory_processor* pr){
    proc = pr;
    proc->add_consumer(this); // Add this consumer to trajectory processor        
}


void Consumer_base::run_in_thread(boost::shared_ptr<Message_channel<boost::shared_ptr<Data_container> > > chan){
    // Call user pre-process
    pre_process_handler();

    boost::shared_ptr<Data_container> data;

    while(chan->recieve(data)){
        consume_frame(data);
    }

    // If we are here than dispatcher thread sent a stop to the queue
    // Consume all remaining frames
    while(!chan->empty()){
        chan->recieve(data);
        consume_frame(data);
    }

    // Call user post-process
    post_process_handler(data->frame_info);
}

void Consumer_base::consume_frame(boost::shared_ptr<Data_container> &data){
    process_window_info(data->frame_info);
    process_frame_data(data->frame);
    process_frame_handler(data->frame_info);
}

void Consumer_base::pre_process(){
}

void Consumer_base::process_frame(const Frame_info& info){

}

void Consumer_base::post_process(const Frame_info& info){
}

void Consumer_base::window_started(const Frame_info& info){
}


void Consumer_base::window_finished(const Frame_info& info){
}


void Consumer_base::process_window_info(Frame_info& info){
    if(info.valid_frame==0){
        // This is the start of the very first window
        win_num = 0; // Set counter
        win_start_time = info.absolute_time;;
        win_start_frame = info.valid_frame;;
        info.win_num = 0;
        info.win_start_frame = info.valid_frame;
        info.win_start_time = info.absolute_time;
        info.win_last_frame = info.valid_frame;
        info.win_last_time = info.absolute_time;
        // Call function for processing
        window_started_handler(info);
    } else {
        // Check the end of window
        if(
            (info.win_size_frames>=0 && info.valid_frame-win_start_frame>=info.win_size_frames)
            ||
            (info.win_size_time>=0 && info.absolute_time-win_start_time>=info.win_size_time)
        ){                        
            info.win_last_frame = info.valid_frame;
            info.win_last_time = info.absolute_time;

            // Previous window finished! Call function for processing.
            // It will see old win_start_time and win_num as needed.
            info.win_num = win_num;
            info.win_start_frame = win_start_frame;
            info.win_start_time = win_start_time;
            window_finished_handler(info);
            // Start new window
            ++win_num;
            // Update window start time
            win_start_frame = info.valid_frame;
            win_start_time = info.absolute_time;

            // Call function for processing
            info.win_num = win_num;
            info.win_start_frame = win_start_frame;
            info.win_start_time = win_start_time;
            window_started_handler(info);
        } else {
            // Window if not finished yet. Just update last_time
            info.win_num = win_num;
            info.win_start_frame = win_start_frame;
            info.win_start_time = win_start_time;
            info.win_last_frame = info.valid_frame;
            info.win_last_time = info.absolute_time;            
        }
    }
}
