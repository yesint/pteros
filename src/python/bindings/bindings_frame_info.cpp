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

#include "bindings_frame_info.h"
#include "pteros/python/bindings_util.h"

using namespace boost::python;
using namespace pteros;

void _set(Frame_info* f, boost::python::tuple data){
    f->absolute_frame = extract<int>(data[0]);
    f->absolute_time = extract<float>(data[1]);
    f->valid_frame = extract<int>(data[2]);
    f->first_time = extract<float>(data[3]);
    f->last_time = extract<float>(data[4]);
    f->first_frame = extract<int>(data[5]);
    f->last_frame = extract<int>(data[6]);
}

void make_bindings_Frame_info(){

    class_<Frame_info>("Frame_info", init<>())
        .def_readonly("absolute_frame",&Frame_info::absolute_frame)
        .def_readonly("absolute_time",&Frame_info::absolute_time)
        .def_readonly("first_frame",&Frame_info::first_frame)
        .def_readonly("first_time",&Frame_info::first_time)
        .def_readonly("last_frame",&Frame_info::last_frame)
        .def_readonly("last_time",&Frame_info::last_time)
        .def_readonly("valid_frame",&Frame_info::valid_frame)
        .def("_set",&_set)
    ;
}
