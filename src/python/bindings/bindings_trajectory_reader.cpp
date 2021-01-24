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
        .def("set_pbc",&Jump_remover::set_pbc)
        .def("set_unwrap_dist",&Jump_remover::set_unwrap_dist)
        .def("set_pbc_atom",&Jump_remover::set_pbc_atom)
    ;

}




