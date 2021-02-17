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


#include "pteros/analysis/task_base.h"
#include <thread>
#include "message_channel.h"
#include "pteros/core/pteros_error.h"
#include "data_container.h"
#include <iostream>

namespace pteros {

using Data_channel = Message_channel<std::shared_ptr<pteros::Data_container> > ;
using Data_channel_ptr = std::shared_ptr<Data_channel> ;

class Task_driver {
public:
    Task_driver(Task_base* _task);
    virtual ~Task_driver();
    void set_data_channel_and_system(const Data_channel_ptr& ch, const System &sys);
    void process_until_end();
    void process_until_end_in_thread ();
    void join_thread();
private:
    Data_channel_ptr channel;
    Task_base* task;
    std::shared_ptr<Data_container> data;
    std::thread t;    
    bool stop_now; // Emergency stop flag for thread
    bool pre_process_done;
};


}
