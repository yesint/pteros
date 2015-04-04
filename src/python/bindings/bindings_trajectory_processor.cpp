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

#include "bindings_trajectory_processor.h"
#include "trajectory_processor_wrapper.h"
#include "consumer_wrapper.h"
#include "pteros/python/bindings_util.h"
#include "pteros/analysis/jump_remover.h"
#include <Eigen/Core>

using namespace pteros;
using namespace Eigen;
using namespace boost::python;
/*
 Bindings to Trajectory_processor could not be direct because
 Python has GIL and can't run the code in several threads as C++ does.
 Instead we make wrapper over Trajectory_processor, which contains exactly one
 consumer inside and exposes same pre_process, process_frame and post_process as Consumer
 This reads file asynchronously, but any other paralellization should be done in python.
 */

// Utility class, which allows overriding pre_process, process_frame and post_process in Python
// and calling that overriden code from C++

/* Here we are forced to work with boost::shared_ptr.
 * We can't use std::shared_ptr since boost.python doesn't wrap it correctly for now...
 */

class _callback: public Trajectory_processor_wrapper {
public:
    _callback(PyObject* p) {
        self = p;
    }

    _callback(PyObject* p, const Options& opt): Trajectory_processor_wrapper(opt) {
        self = p;
    }

    void pre_process() {        
        call_method<void>(self, "pre_process");
    }

    void process_frame(const Frame_info& info) {
        call_method<void>(self, "process_frame", info);
    }

    void post_process(const Frame_info& info) {
        call_method<void>(self, "post_process", info);
    }

    PyObject* self;

};

void jump_remover_dims(Jump_remover* r, PyObject* dims_to_wrap){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims_to_wrap)
    r->set_dimensions(dim);
}


void make_bindings_Trajectory_processor(){

    // Need to wrap base class first
    class_<Trajectory_processor, boost::noncopyable>("Trajectory_processor_base", init<>())
    ;

    class_<Trajectory_processor_wrapper, boost::noncopyable, boost::shared_ptr<_callback> , bases<Trajectory_processor> >("Trajectory_processor", init<>())
        .def(init<const Options&>() )
        .def("set_options",&Trajectory_processor_wrapper::set_options)
        .def("run",&Trajectory_processor_wrapper::run)
        .def("get_system",&Trajectory_processor_wrapper::get_system,return_value_policy<reference_existing_object>())
        .def("get_frame_ptr",&Trajectory_processor_wrapper::get_frame_ptr,return_value_policy<reference_existing_object>())
        .def("initialize",&Trajectory_processor_wrapper::initialize)
        .def("help",&Trajectory_processor_wrapper::help)
    ;

    // Also wrap Jump_remover. It will be created for each python plugin on python side
    class_<Jump_remover>("Jump_remover",init<>())
        .def("add_atoms",&Jump_remover::add_atoms)
        .def("set_dimensions",&jump_remover_dims)
        .def("set_unwrap_dist",&Jump_remover::set_unwrap_dist)
        .def("remove_jumps",&Jump_remover::remove_jumps)
    ;

}
