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

#ifndef TASK_PYTHON_H
#define TASK_PYTHON_H

#include "pteros/analysis/selection_analysis.h"

#include <boost/python.hpp>
#include <numpy/noprefix.h>
#include <boost/python/numeric.hpp>


namespace pteros {

class Task_python: public Consumer {
    //friend class Proxy;
public:
    Task_python(Trajectory_processor* engine, Options_tree* opt);

    virtual void pre_process();
    virtual bool process_frame(const Frame_info& info);
    virtual void post_process(const Frame_info& info);
    //virtual void window_started_slot(const Trajectory_processor_info& info);
    //virtual void window_finished_slot(const Trajectory_processor_info& info);
    static void print_help();

private:
    Options_tree* options;
    boost::python::object main_namespace;   

    template<class T>
    void fetch_infile_option(string option){
        using namespace boost::python;
        try {
            options->add_value("trajectory/"+option,
                               extract<T>(main_namespace["trajectory"].attr(option.c_str()))());
            T val = extract<T>(main_namespace["trajectory"].attr(option.c_str()))();
            if(val!=-1)
                cout << "trajectory."+option+" = " << val << endl;
        } catch(error_already_set const &) {
            PyErr_Clear();
        }
    }

};

}

#endif
