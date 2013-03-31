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

#include "consumer_wrapper.h"
#include "trajectory_processor_wrapper.h"

using namespace pteros;
using namespace boost::python;

void Consumer_wrapper::pre_process(){
    try {        
        dynamic_cast<Trajectory_processor_wrapper*>(proc)->pre_process();
    } catch (error_already_set& e){
        PyErr_Print();
    }
}

bool Consumer_wrapper::process_frame(const Frame_info& info){    
    try {
        return dynamic_cast<Trajectory_processor_wrapper*>(proc)->process_frame(info);
    } catch (error_already_set& e){
        PyErr_Print();
        return false;
    }
}

void Consumer_wrapper::post_process(const Frame_info& info){
    try {        
        dynamic_cast<Trajectory_processor_wrapper*>(proc)->post_process(info);
    } catch (error_already_set& e){
        PyErr_Print();
    }
}

