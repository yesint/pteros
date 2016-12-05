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

struct Frame_info_suite : boost::python::pickle_suite
{
    static boost::python::tuple getstate(boost::python::object w_obj)
    {
        Frame_info const& w = boost::python::extract<Frame_info const&>(w_obj)();
        return boost::python::make_tuple(
                  w_obj.attr("__dict__"), //If the python object has other attributes, they will be stored in the dict
                  w.absolute_frame,
                  w.absolute_time,
                  w.valid_frame,
                  w.first_time,
                  w.last_time,
                  w.first_frame,
                  w.last_frame);
    }

    static void setstate(boost::python::object w_obj, boost::python::tuple state)
    {
        using namespace boost::python;
        Frame_info& w = extract<Frame_info&>(w_obj)();
        // restore the object's __dict__
        dict d = extract<dict>(w_obj.attr("__dict__"))();
        d.update(state[0]);
        w.absolute_frame = extract<int>(state[1]);
        w.absolute_time = extract<float>(state[2]);
        w.valid_frame = extract<int>(state[3]);
        w.first_time = extract<float>(state[4]);
        w.last_time = extract<float>(state[5]);
        w.first_frame = extract<int>(state[6]);
        w.last_frame = extract<int>(state[7]);
    }
    static bool getstate_manages_dict() { return true; }
};

void make_bindings_Frame_info(){

    class_<Frame_info>("Frame_info", init<>())
        .def_readonly("absolute_frame",&Frame_info::absolute_frame)
        .def_readonly("absolute_time",&Frame_info::absolute_time)
        .def_readonly("first_frame",&Frame_info::first_frame)
        .def_readonly("first_time",&Frame_info::first_time)
        .def_readonly("last_frame",&Frame_info::last_frame)
        .def_readonly("last_time",&Frame_info::last_time)
        .def_readonly("valid_frame",&Frame_info::valid_frame)
        .enable_pickling()
        .def_pickle(Frame_info_suite())
    ;
}
