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

#ifndef COMPILED_PLUGIN_H
#define COMPILED_PLUGIN_H

#include "pteros/analysis/task_plugin.h"

#ifndef STANDALONE_PLUGINS

// Make a Python extension module from this plugin

#include "pteros/python/bindings_util.h"

using namespace boost::python;

#define CREATE_COMPILED_PLUGIN(_name,_dum) \
void Pteros_error_translator(const pteros::Pteros_error& e) { \
  PyErr_SetString(PyExc_UserWarning, const_cast<pteros::Pteros_error&>(e).what()); \
} \
_name* get_instance_##_name(const Options& opt){ return new _name(opt); }\
BOOST_PYTHON_MODULE(_name) \
{ \
    import_array(); \
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray"); \
    register_exception_translator<pteros::Pteros_error>(&Pteros_error_translator);\
    class_<_name,bases<Task_plugin>,boost::noncopyable>("Task", init<const pteros::Options&>())\
    .def("help",&_name::help) \
    ; \
    def("get_instance",&get_instance_##_name, return_value_policy<reference_existing_object>());\
}

#else //STANDALONE_PLUGINS

// Make a stand-alone executable from this plugin

#include "pteros/analysis/trajectory_reader.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

using namespace pteros;
using namespace std;

// Skeleton of the driver program

#define CREATE_COMPILED_PLUGIN(_name, _collector)\
int main(int argc, char** argv){\
    try {\
        Options options;\
        parse_command_line(argc,argv,options);\
        Trajectory_reader engine(options);\
        auto task = new _name(options);\
        engine.add_task(task);\
        engine.register_collector(_collector);\
        cout << "-------------------------------------------------------------" << endl;\
        cout << "  This is stand-alone Pteros analysis plugin '" #_name "'" << endl;\
        cout << "-------------------------------------------------------------" << endl;\
        if(!options.has("f") && !options.has("help")){\
            cout << "Usage:" << endl;\
            cout << "\tpteros_" #_name " -f <files> <task options>" << endl;\
            cout << "\n\tFor specific task options use '-help task'" << endl;\
            cout << "\tFor trajectory processing options use '-help traj'" << endl;\
            cout << "\tFor all available options use '-help all' or just '-help'" << endl;\
            return 1;\
        }\
        if(options.has("help")){\
            string help = options("help","").as_string();\
            if(help=="traj"){\
                cout << engine.help() << endl;\
            } else if(help=="task"){\
                cout << task->help() << endl;\
            } else {\
                cout << task->help() << endl << endl;\
                cout << engine.help() << endl;\
            }\
            return 1;\
        }\
        engine.run();\
    } catch (const std::exception& e) { \
        LOG()->error(e.what()); \
    } catch(...) { \
        LOG()->error("Unknown error"); \
    }\
}


#endif //STANDALONE_PLUGINS

#endif //COMPILED_PLUGIN_H
