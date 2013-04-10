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

#include "trajectory_processor_wrapper.h"

using namespace pteros;

void Trajectory_processor_wrapper::initialize(){
    cons_p = boost::shared_ptr<Consumer_wrapper>(new Consumer_wrapper(this));    
}

Trajectory_processor_wrapper::Trajectory_processor_wrapper(): Trajectory_processor(){
    initialize();
}
Trajectory_processor_wrapper::Trajectory_processor_wrapper(Options_tree& opt): Trajectory_processor(opt){
    initialize();
}

System* Trajectory_processor_wrapper::get_system(){
    return cons_p->get_system();
}

Frame *Trajectory_processor_wrapper::get_frame_ptr(){
    return cons_p->get_frame_ptr();
}
