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

#include "bindings_trajectory_processor.h"
#include "trajectory_processor_wrapper.h"
#include "consumer_wrapper.h"

/*
 Bindings to Trajectory_processor could not be direct because
 Python has GIL and can't run the code in several threads as C++ does.
 Instead we make wrapper over Trajectory_processor, which contains exactly one
 consumer inside and exposes same pre_process, process_frame and post_process as Consumer
 This reads file asynchronously, but any other paralellization should be done in python.
 */

// Utility class, which allows overriding pre_process, process_frame and post_process in Python
// and calling that overriden code from C++
class _callback: public Trajectory_processor_wrapper {
public:
    _callback(PyObject* p) {
        self = p;
    }

    _callback(PyObject* p, Options_tree& opt): Trajectory_processor_wrapper(opt) {
        self = p;
    }

    void pre_process() {        
        call_method<void>(self, "pre_process");
    }

    bool process_frame(const Frame_info& info) {
        return call_method<bool>(self, "process_frame", info);        
    }

    void post_process(const Frame_info& info) {
        call_method<void>(self, "post_process", info);
    }

    PyObject* self;
};


void make_bindings_Trajectory_processor(){

    class_<Trajectory_processor_wrapper, boost::noncopyable, boost::shared_ptr<_callback> >("Trajectory_processor", init<>())
        .def(init<Options_tree&>() )
        .def("set_options",&Trajectory_processor_wrapper::set_options)
        .def("run",&Trajectory_processor_wrapper::run)
        .def("get_system",&Trajectory_processor_wrapper::get_system,return_value_policy<reference_existing_object>())
    ;

    // Exposed only to be created in the script and bound to trajectory processor
    // Can't inherit or execute any Python code
    class_<Consumer>("Compiled_task",init<Trajectory_processor*>())
    ;
}
