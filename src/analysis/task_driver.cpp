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


#include "task_driver.h"

using namespace std;
using namespace pteros;

TaskDriver::TaskDriver(TaskBase *_task): task(_task), stop_now(false)
{
    //cout << "ctor: Task_driver" << endl;
}

TaskDriver::~TaskDriver()
{
    if(t.joinable()){
        // Stop the thread
        stop_now = true;
        task->log->error("Ups! Stopping task driver thread on outer exception...");
        t.join();
    }
}

void TaskDriver::set_data_channel_and_system(const DataChannel_ptr &ch, const System &sys){
    channel = ch;
    task->put_system(sys);
}

void TaskDriver::process_until_end() {
    pre_process_done = false;
    while(channel->recieve(data)){
        if(stop_now) return; // Emergency stop point

        task->put_frame(data->frame);
        if(!pre_process_done){
            task->pre_process_handler();
            pre_process_done = true;
        }
        task->process_frame_handler(data->frame_info);
        ++task->n_consumed;
    }
    if(task->n_consumed>0){        
        task->post_process_handler(data->frame_info);        
    } else {
        task->log->warn("No frames consumed!");
    }
}

void TaskDriver::process_until_end_in_thread() {
    t = std::thread(&TaskDriver::process_until_end, this);
}

void TaskDriver::join_thread() { t.join(); }

