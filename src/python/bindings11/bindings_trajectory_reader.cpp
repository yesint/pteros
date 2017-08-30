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

#include "pteros/analysis/trajectory_reader.h"
#include "pteros/analysis/options.h"
#include "pteros/analysis/task_plugin.h"
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>


namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace pybind11::literals;

class Task_base_py : public Task_plugin {
public:
    /* Inherit the constructors */
    Task_base_py(const Options& opt): Task_plugin(opt) { }
    virtual void pre_process_pub() = 0;
    virtual void process_frame_pub(const Frame_info& info) = 0;
    virtual void post_process_pub(const Frame_info& info) = 0;


protected:
    virtual bool is_parallel(){return false;}
    // Call corresponding public functions
    virtual void pre_process(){ pre_process_pub(); }
    virtual void process_frame(const Frame_info& info){ process_frame_pub(info); }
    virtual void post_process(const Frame_info& info){ post_process_pub(info); }

};

class Task_py : public Task_base_py {
public:
    /* Inherit the constructors */
    using Task_base_py::Task_base_py;

    /* Trampoline (need one for each virtual function) */
    void pre_process_pub() override {
        PYBIND11_OVERLOAD_PURE(
            void, /* Return type */
            Task_base_py,      /* Parent class */
            pre_process_pub          /* Name of function in C++ (must match Python name) */
                          /* Argument(s) */
        );
    }

    void process_frame_pub(const Frame_info& info) override {
        PYBIND11_OVERLOAD_PURE(
            void, /* Return type */
            Task_base_py,      /* Parent class */
            process_frame_pub,          /* Name of function in C++ (must match Python name) */
            info              /* Argument(s) */
        );
    }

    void post_process_pub(const Frame_info& info) override {
        PYBIND11_OVERLOAD_PURE(
            void, /* Return type */
            Task_base_py,      /* Parent class */
            post_process_pub,          /* Name of function in C++ (must match Python name) */
            info              /* Argument(s) */
        );
    }

};


void make_bindings_Trajectory_reader(py::module& m){

    py::class_<Task_base_py,Task_py>(m, "Task_base")
        .def(py::init<const Options&>())
        .def("pre_process",&Task_base_py::pre_process_pub)
        .def("process_frame",&Task_base_py::process_frame_pub)
        .def("post_process",&Task_base_py::post_process_pub)
    ;

    py::class_<Trajectory_reader>(m, "Trajectory_reader")
        .def(py::init<>())
        .def(py::init<const Options&>())
        .def("set_options",&Trajectory_reader::set_options)
        .def("help",&Trajectory_reader::help)
        .def("run",&Trajectory_reader::run)
        .def("add_task",[](Trajectory_reader* r, const py::object& o){
                auto p = o.cast<Task_base_py*>();
                r->add_task(p);
            })
    ;

}
