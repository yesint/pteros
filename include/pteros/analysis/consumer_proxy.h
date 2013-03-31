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


#ifndef CONSUMER_PROXY_H
#define CONSUMER_PROXY_H

#include "pteros/analysis/consumer_base.h"

namespace pteros {

/** Proxy consumer for python bindings
  */
class Consumer_proxy: public Consumer_base {
public:
    Consumer_proxy(Trajectory_processor* pr): Consumer_base(pr){
        proxy_frame_ptr = NULL;
    }

    Frame* get_frame_ptr(){
        return proxy_frame_ptr;
    }

protected:
    Frame* proxy_frame_ptr;

    virtual void process_frame_data(Frame& data);
};

}

#endif
