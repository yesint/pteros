/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
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



class Task_py : public TaskPlugin {
public:
    /* Inherit the constructors */    
    using TaskPlugin::TaskPlugin;

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
            TaskPlugin,      /* Parent class */
            pre_process,          /* Name of function in C++ (must match Python name) */
                          /* Argument(s) */
        );
    }

    void process_frame(const FrameInfo& info) override {
        py::gil_scoped_acquire acquire;
        PYBIND11_OVERLOAD_PURE(
            void, /* Return type */
            TaskPlugin,      /* Parent class */
            process_frame,          /* Name of function in C++ (must match Python name) */
            info              /* Argument(s) */
        );
    }

    void post_process(const FrameInfo& info) override {
        py::gil_scoped_acquire acquire;
        PYBIND11_OVERLOAD_PURE(
            void, /* Return type */
            TaskPlugin,      /* Parent class */
            post_process,          /* Name of function in C++ (must match Python name) */
            info              /* Argument(s) */
        );
    }
protected:
    bool is_parallel() override { return false; }
};


void make_bindings_Trajectory_reader(py::module& m){

    py::class_<TaskPlugin,Task_py,std::shared_ptr<TaskPlugin>>(m, "TaskBase")
        .def(py::init<const Options&>())        

        .def("pre_process",&TaskPlugin::pre_process)
        .def("process_frame",&TaskPlugin::process_frame)
        .def("post_process",&TaskPlugin::post_process)

        .def_readonly("system",&TaskPlugin::system)
        .def_property_readonly("id",&TaskPlugin::get_id)
        .def_readonly("jump_remover",&TaskPlugin::jump_remover)
        .def_readonly("options",&TaskPlugin::options)
        .def_property_readonly("log",[](TaskPlugin* obj){return obj->log.get();},py::return_value_policy::reference_internal)

        .def_property("_class_name",[](Task_py* obj){return obj->_class_name;}, [](Task_py* obj, const string& s){obj->_class_name=s;})
    ;

    py::class_<TrajectoryReader>(m, "TrajectoryReader")
        .def(py::init<>())
        .def(py::init<const Options&>())
        .def("set_options",&TrajectoryReader::set_options)
        .def("help",&TrajectoryReader::help)
        // We release GIL before starting run and will acquire it in individual tasks
        .def("run",[](TrajectoryReader* r){
            py::gil_scoped_release release;
            r->run();
        })

        .def("add_task",[](TrajectoryReader* r, const py::object& o){
            auto p = o.cast<shared_ptr<TaskPlugin>>();
            r->add_task(p);
        })
    ;

    py::class_<FrameInfo>(m,"FrameInfo")
        .def_readonly("absolute_frame",&FrameInfo::absolute_frame)
        .def_readonly("absolute_time",&FrameInfo::absolute_time)
        .def_readonly("first_frame",&FrameInfo::first_frame)
        .def_readonly("first_time",&FrameInfo::first_time)
        .def_readonly("last_frame",&FrameInfo::last_frame)
        .def_readonly("last_time",&FrameInfo::last_time)
        .def_readonly("valid_frame",&FrameInfo::valid_frame)
    ;

    py::class_<JumpRemover>(m,"JumpRemover")
        .def("add_atoms",&JumpRemover::add_atoms)
        .def("set_pbc",&JumpRemover::set_pbc)
        .def("set_unwrap_dist",&JumpRemover::set_unwrap_dist)
        .def("set_pbc_atom",&JumpRemover::set_pbc_atom)
    ;

}




