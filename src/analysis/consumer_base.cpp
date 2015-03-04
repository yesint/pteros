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
#include "pteros/core/pteros_error.h"

using namespace std;
using namespace pteros;

Consumer_base::Consumer_base(Trajectory_processor* pr){
    proc = pr;
    proc->add_consumer(this); // Add this consumer to trajectory processor
}

void Consumer_base::run_in_thread(std::shared_ptr<Message_channel<std::shared_ptr<pteros::Data_container> > > &chan){
    // pre-process is called on very first valid frame in consume_frame()

    std::shared_ptr<Data_container> data;

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

void Consumer_base::consume_frame(std::shared_ptr<Data_container> &data){
    // Check number of atoms
    if(data->frame.coord.size()!=system.num_atoms()){
        throw Pteros_error("Wrong number of atoms in the trajectory frame!");
    }
    // If this is the very first valid frame call pre_process
    if(data->frame_info.valid_frame==0){
        pre_process_handler();
    }

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

