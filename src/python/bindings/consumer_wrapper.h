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

#ifndef CONSUMER_WRAPPER_H
#define CONSUMER_WRAPPER_H

#include "pteros/python/bindings_util.h"
#include "consumer_proxy.h"

namespace pteros {

class Consumer_wrapper: public Consumer_proxy {
public:
    Consumer_wrapper(Trajectory_processor* pr): Consumer_proxy(pr){
    }

    ~Consumer_wrapper(){
    }

protected:
    virtual void pre_process();
    virtual void process_frame(const Frame_info& info);
    virtual void post_process(const Frame_info &info);    
};

}
#endif
