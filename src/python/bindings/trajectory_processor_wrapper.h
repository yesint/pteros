/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2014, Semen Yesylevskyy
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

#ifndef TRAJECTORY_PROCESSOR_WRAPPER_H
#define TRAJECTORY_PROCESSOR_WRAPPER_H

#include "pteros/python/bindings_util.h"
#include "pteros/analysis/trajectory_processor.h"
#include "consumer_wrapper.h"

namespace pteros {

class Trajectory_processor_wrapper: public Trajectory_processor {
public:
    void initialize();

    Trajectory_processor_wrapper();
    Trajectory_processor_wrapper(Options_tree& opt);
    ~Trajectory_processor_wrapper(){}    

    virtual void pre_process() = 0;
    virtual void process_frame(const Frame_info& info) = 0;
    virtual void post_process(const Frame_info& info) = 0;

    System* get_system();   
    // Gets pointer to internal frame obtained from reader
    Frame* get_frame_ptr();

private:
    boost::shared_ptr<Consumer_wrapper> cons_p;
};

}

#endif
