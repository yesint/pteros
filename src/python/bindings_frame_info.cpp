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

#include "bindings_frame_info.h"

void make_bindings_Frame_info(){

    class_<Frame_info>("Frame_info", init<>())
        .def_readwrite("absolute_frame",&Frame_info::absolute_frame)
        .def_readwrite("first_frame",&Frame_info::first_frame)
        .def_readwrite("first_time",&Frame_info::first_time)
        .def_readwrite("last_frame",&Frame_info::last_frame)
        .def_readwrite("last_time",&Frame_info::last_time)
        .def_readwrite("valid_frame",&Frame_info::valid_frame)

        .def_readwrite("win_size_frames",&Frame_info::win_size_frames)
        .def_readwrite("win_size_time",&Frame_info::win_size_time)
        .def_readwrite("win_start_time",&Frame_info::win_start_time)
        .def_readwrite("win_start_frame",&Frame_info::win_start_frame)
        .def_readwrite("win_last_time",&Frame_info::win_last_time)
        .def_readwrite("win_last_frame",&Frame_info::win_last_frame)
        .def_readwrite("win_num",&Frame_info::win_num)
    ;
}
