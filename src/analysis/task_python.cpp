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

#include "task_python.h"
#include <fstream>

using namespace std;
using namespace pteros;
using namespace boost;
using namespace boost::python;
using namespace boost::python::numeric;

Task_python::Task_python(Trajectory_processor *engine, Options_tree *opt): Consumer(engine){
    options = opt;
    try {
        cout << "Initializing embedded python interpreter..." << endl;
        Py_Initialize();
        object main_module = import("__main__");
        main_namespace = main_module.attr("__dict__");

        // Import a pteros_py module
        // EXTRA_PYTHON_PATH macro is defined in CMake and pathed by -D... flag
#ifdef EXTRA_PYTHON_PATH
        string module_path = options->get_value<string>("module_path",EXTRA_PYTHON_PATH());
#else
        string module_path = options->get_value<string>("module_path","");
#endif
        cout << "Extra module path is " << module_path << endl;
        string todo = "import sys\n"
                "sys.path.append('"+module_path+"')\n"
                "print 'Importing pteros_py module into embedded interpreter...'\n"
                "from pteros_py import *\n"
                "print 'Checking import - the System class is imported as:',System\n"
                ;
        exec(todo.c_str(),main_namespace);

        main_namespace["system"] = ptr(&system);
        main_namespace["id"] = id;

        todo =  "class Trajectory:\n"
                "   def __init__(self):\n"
                "       self.files = []\n"
                "       self.first_frame = -1\n"
                "       self.last_frame = -1\n"
                "\n"
                "trajectory = Trajectory()\n"
                ;
        exec(todo.c_str(),main_namespace);

        string script_file = options->get_value<string>("script");
        cout << "Loading user-defined script '"+ script_file +"'..." << endl;
        exec_file(script_file.c_str(),main_namespace);

        cout << "Processing in-script trajectory options (if any)..." << endl;

        try {
            boost::python::list l = extract<boost::python::list>(main_namespace["trajectory"]);
            for(int i=0; i<boost::python::len(l); ++i){
                options->add_value("trajectory",extract<string>(l[i])());
                cout << "trajectory += " << extract<string>(l[i])() << endl;
            }
        } catch(error_already_set const &) {
            cout << "No in-script trajectory options" << endl;
            PyErr_Clear();
        }

        try {
            options->add_value("trajectory/first_frame",
                               extract<int>(main_namespace["trajectory"].attr("first_frame"))());
            cout << "trajectory.first_frame = "
                 << extract<int>(main_namespace["trajectory"].attr("first_frame"))() << endl;
        } catch(error_already_set const &) {
            cout << "No first_frame, well" << endl;
            PyErr_Clear();
        }

        // Create a task instance
        exec("task = Task()",main_namespace);


    } catch(error_already_set const &) {
        PyErr_Print();
        throw Pteros_error("Failed to initialize embedded Python interpreter!");
    }
}

void Task_python::pre_process(){
    try {
        exec("task.pre_process()",main_namespace);
    } catch(error_already_set const &) {
         PyErr_Print();
         throw Pteros_error("Python pre_process failed!");
    }
}

bool Task_python::process_frame(const Frame_info &info){
    bool ok;
    try {
        main_namespace["_info"] = ptr(&info);
        exec("_ret = task.process_frame(_info)",main_namespace);
        ok = extract<bool>(main_namespace["_ret"]);
    } catch(error_already_set const &) {
         PyErr_Print();
         throw Pteros_error("Python process_frame failed!");
    }
    return ok;
}

void Task_python::post_process(const Frame_info &info){
    try {        
        main_namespace["_info"] = ptr(&info);
        exec("task.post_process(_info)" , main_namespace);

    } catch(error_already_set const &) {
         PyErr_Print();
         throw Pteros_error("Python post_process failed!");
    }
}

//void Task_box::window_started_slot(const Trajectory_processor_info& info){}
//void Task_box::window_finished_slot(const Trajectory_processor_info& info){}

void Task_python::print_help(){
    cout << "Task python:\n"
            "-----------------\n"
            "Executes user python script.\n"
         << endl;
}
