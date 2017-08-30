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

#include "bindings_trajectory_reader.h"
#include "pteros/python/bindings_util.h"
#include "pteros/analysis/trajectory_reader.h"
#include "pteros/analysis/task_plugin.h"
#include "pteros/core/logging.h"
#include <Eigen/Core>

using namespace pteros;
using namespace Eigen;
using namespace boost::python;
using namespace std;


class Trajectory_reader_adaptor {
public:

    TASK_SERIAL(Task_inner)
    public:
        Task_inner(Trajectory_reader_adaptor* r): reader(r), Task_base() {}
    protected:
        virtual void pre_process() {
            try {
                reader->pre_process_cb(&system);
            } catch (error_already_set& e){
                LOG()->error("Error occured in Python code:");
                PyErr_Print();
                exit(1);
            }
        }

        virtual void process_frame(const Frame_info& info) {
            try {
               reader->process_frame_cb(info);
            } catch (error_already_set& e){
                LOG()->error("Error occured in Python code:");
                PyErr_Print();
                exit(1);
            }
        }

        virtual void post_process(const Frame_info& info) {
            try {
                reader->post_process_cb(info);
            } catch (error_already_set& e){
                LOG()->error("Error occured in Python code:");
                PyErr_Print();
                exit(1);
            }
        }
    private:
        Trajectory_reader_adaptor* reader;
    };

    Trajectory_reader_adaptor(){

    }

    Trajectory_reader_adaptor(const Options& opt,
                              boost::python::object pre_process_handler,
                              boost::python::object process_frame_handler,
                              boost::python::object post_process_handler) {
        reader.set_options(opt);

        pre_process_cb = [pre_process_handler](System* sys) { pre_process_handler(ptr(sys)); };
        process_frame_cb = [process_frame_handler](const Frame_info& info) { process_frame_handler(info); };
        post_process_cb = [post_process_handler](const Frame_info& info) { post_process_handler(info); };
    }

    void add_python_tasks_runner(){
        reader.add_task( new Task_inner(this) );
    }

    void run(){
        reader.run();
    }

    void add_task(object& obj){
        Task_plugin* ptr = extract<Task_plugin*>(obj);
        reader.add_task( dynamic_cast<Task_base*>(ptr) );
    }

    std::function<void(System*)> pre_process_cb;
    std::function<void(const Frame_info&)> process_frame_cb;
    std::function<void(const Frame_info&)> post_process_cb;

private:
    Trajectory_reader reader;
};

class Logger {
public:
    Logger(const string& name){
        log = std::make_shared<spdlog::logger>(name, Log::instance().console_sink);
        log->set_pattern(Log::instance().generic_pattern);
    }

    void info(const string& msg){ log->info(msg); }
    void error(const string& msg){ log->error(msg); }
    void warn(const string& msg){ log->warn(msg); }
    void set_pattern(const string& pat){ log->set_pattern(pat); }

private:
    std::shared_ptr<spdlog::logger> log;
};


void jump_remover_dims(Jump_remover* r, PyObject* dims_to_wrap){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims_to_wrap)
    r->set_dimensions(dim);
}



void make_bindings_Trajectory_reader(){

    class_<Task_plugin, boost::noncopyable>("Task_plugin", no_init)
    ;

    class_<Trajectory_reader_adaptor, boost::noncopyable>("Trajectory_reader", init<>())
        .def(init<const Options&,boost::python::object,boost::python::object,boost::python::object >() )        
        .def("run", &Trajectory_reader_adaptor::run)
        .def("add_task", &Trajectory_reader_adaptor::add_task)
        .def("add_python_tasks_runner", &Trajectory_reader_adaptor::add_python_tasks_runner)
    ;

    // Bindings for logger
    class_<Logger>("Logger", init<const string&>())
        .def("info",&Logger::info)
        .def("warn",&Logger::warn)
        .def("error",&Logger::error)
        .def("set_pattern",&Logger::set_pattern)
    ;

    // Also wrap Jump_remover. It will be created for each python plugin on python side
    class_<Jump_remover>("Jump_remover",init<>())
        .def("add_atoms",&Jump_remover::add_atoms)
        .def("set_dimensions",&jump_remover_dims)
        .def("set_unwrap_dist",&Jump_remover::set_unwrap_dist)
        .def("remove_jumps",&Jump_remover::remove_jumps)
        .def("set_leading_index",&Jump_remover::set_leading_index)
    ;


}
