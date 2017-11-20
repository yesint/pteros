/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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
#include "bindings_util.h"

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace pybind11::literals;



class Task_py : public Task_plugin {
public:
    /* Inherit the constructors */    
    using Task_plugin::Task_plugin;

    // Workaround to set the correct class name for logger
    // This is filled on Python side by the name of derived class
    string _class_name;

    void set_id(int id) override {                
        log = create_logger(fmt::format("{}.{}",_class_name,id));
        task_id = id;        
    }

    /* Trampoline (need one for each virtual function) */
    void pre_process() override {
        py::gil_scoped_acquire acquire;
        PYBIND11_OVERLOAD_PURE(
            void, /* Return type */
            Task_plugin,      /* Parent class */
            pre_process,          /* Name of function in C++ (must match Python name) */
                          /* Argument(s) */
        );        
    }

    void process_frame(const Frame_info& info) override {
        py::gil_scoped_acquire acquire;
        PYBIND11_OVERLOAD_PURE(
            void, /* Return type */
            Task_plugin,      /* Parent class */
            process_frame,          /* Name of function in C++ (must match Python name) */
            info              /* Argument(s) */
        );        
    }

    void post_process(const Frame_info& info) override {
        py::gil_scoped_acquire acquire;
        PYBIND11_OVERLOAD_PURE(
            void, /* Return type */
            Task_plugin,      /* Parent class */
            post_process,          /* Name of function in C++ (must match Python name) */
            info              /* Argument(s) */
        );
    }
protected:
    bool is_parallel() override { return false; }
};


void make_bindings_Trajectory_reader(py::module& m){

    py::class_<Task_plugin,Task_py,std::shared_ptr<Task_plugin>>(m, "Task_base")
        .def(py::init<const Options&>())        

        .def("pre_process",&Task_plugin::pre_process)
        .def("process_frame",&Task_plugin::process_frame)
        .def("post_process",&Task_plugin::post_process)

        .def_readonly("system",&Task_plugin::system)
        .def_property_readonly("id",&Task_plugin::get_id)
        .def_readonly("jump_remover",&Task_plugin::jump_remover)
        .def_readonly("options",&Task_plugin::options)
        .def_property_readonly("log",[](Task_plugin* obj){return obj->log.get();},py::return_value_policy::reference_internal)

        .def_property("_class_name",[](Task_py* obj){return obj->_class_name;}, [](Task_py* obj, const string& s){obj->_class_name=s;})
    ;

    py::class_<Trajectory_reader>(m, "Trajectory_reader")
        .def(py::init<>())
        .def(py::init<const Options&>())
        .def("set_options",&Trajectory_reader::set_options)
        .def("help",&Trajectory_reader::help)
        // We release GIL before starting run and will acquire it in individual tasks
        .def("run",[](Trajectory_reader* r){
            py::gil_scoped_release release;
            r->run();
        })

        .def("add_task",[](Trajectory_reader* r, const py::object& o){
            auto p = o.cast<shared_ptr<Task_plugin>>();
            r->add_task(p);
        })
    ;

    py::class_<Frame_info>(m,"Frame_info")
        .def_readonly("absolute_frame",&Frame_info::absolute_frame)
        .def_readonly("absolute_time",&Frame_info::absolute_time)
        .def_readonly("first_frame",&Frame_info::first_frame)
        .def_readonly("first_time",&Frame_info::first_time)
        .def_readonly("last_frame",&Frame_info::last_frame)
        .def_readonly("last_time",&Frame_info::last_time)
        .def_readonly("valid_frame",&Frame_info::valid_frame)
    ;

    py::class_<Jump_remover>(m,"Jump_remover")
        .def("add_atoms",&Jump_remover::add_atoms)
        .def("set_dimensions",&Jump_remover::set_dimensions)
        .def("set_unwrap_dist",&Jump_remover::set_unwrap_dist)
        .def("set_leading_index",&Jump_remover::set_leading_index)
    ;

}
