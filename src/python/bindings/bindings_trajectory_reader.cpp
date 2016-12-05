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

#include "bindings_trajectory_reader.h"
#include "pteros/python/bindings_util.h"
#include "pteros/analysis/trajectory_reader.h"
#include "pteros/analysis/task_base.h"
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
            reader->pre_process_cb(&system);
        }
        virtual void process_frame(const Frame_info& info) {
            reader->process_frame_cb(info);
        }
        virtual void post_process(const Frame_info& info) {
            reader->post_process_cb(info);
        }
    private:
        Trajectory_reader_adaptor* reader;
    };

    Trajectory_reader_adaptor(){}

    Trajectory_reader_adaptor(const Options& opt,
                              boost::python::object pre_process_handler,
                              boost::python::object process_frame_handler,
                              boost::python::object post_process_handler) {
        reader.set_options(opt);

        pre_process_cb = [&pre_process_handler](System* sys) { pre_process_handler(ptr(sys)); };
        process_frame_cb = [&process_frame_handler](const Frame_info& info) { process_frame_handler(info); };
        post_process_cb = [&post_process_handler](const Frame_info& info) { post_process_handler(info); };

        reader.add_task( new Task_inner(this) );
        reader.run();
    }

    std::function<void(System*)> pre_process_cb;
    std::function<void(const Frame_info&)> process_frame_cb;
    std::function<void(const Frame_info&)> post_process_cb;

private:
    Trajectory_reader reader;
};


void make_bindings_Trajectory_reader(){

    class_<Trajectory_reader_adaptor, boost::noncopyable>("Trajectory_reader", init<>())
        .def(init<const Options&,boost::python::object,boost::python::object,boost::python::object >() )
    ;
}
