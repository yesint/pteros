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

#include "pteros/analysis/consumer.h"

using namespace pteros;

void Consumer::process_frame_data(Frame &data){
    system.Frame_data(0) = data;
}

void Consumer::process_frame_handler(const Frame_info &info){
    // Remove jumps
    jump_remover.remove_jumps(system);
    // Call user callback
    process_frame(info);
}

