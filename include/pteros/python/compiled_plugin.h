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

#ifndef STANDALONE_PLUGINS

// Make a Python extension module from this plugin

#include "pteros/python/bindings_util.h"
#include "pteros/python/compiled_plugin_base.h"

using namespace boost::python;

#define CREATE_COMPILED_PLUGIN(_name) \
void Pteros_error_translator(const pteros::Pteros_error& e) { \
  PyErr_SetString(PyExc_UserWarning, const_cast<pteros::Pteros_error&>(e).what().c_str()); \
} \
BOOST_PYTHON_MODULE(_name) \
{ \
    import_array(); \
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray"); \
    register_exception_translator<pteros::Pteros_error>(&Pteros_error_translator); \
    class_<_name,boost::noncopyable>("Task", init<pteros::Trajectory_processor*,pteros::Options_tree*>()) \
    .def_readwrite("label",&_name::label) \
    .def("help",&_name::help) \
    ; \
}

#else //STANDALONE_PLUGINS

// Make a stand-alone executable from this plugin

#include "pteros/core/pteros_error.h"
#include "pteros/python/compiled_plugin_base.h"

using namespace pteros;
using namespace std;

#define CREATE_COMPILED_PLUGIN(_name) \
int main(int argc, char** argv){\
    try {\
        Options_tree options;\
        options.from_command_line(argc,argv);\
        Trajectory_processor engine(options);\
        _name task(&engine,&options);\
        task.label = #_name;\
        cout << "-------------------------------------------------------------" << endl;\
        cout << "  This is stand-alone Pteros analysis plugin '" #_name "'" << endl;\
        cout << "-------------------------------------------------------------" << endl;\
        if(options.count_options("trajectory")==0 && options.count_options("help")==0){\
            cout << "Usage:" << endl;\
            cout << "\tpteros_" #_name " --trajectory[<traj options>] <task options>" << endl;\
            cout << "\n\tFor specific task options use '--help task'" << endl;\
            cout << "\tFor trajectory processing options use '--help traj'" << endl;\
            cout << "\tFor all available options use '--help all' or just '--help'" << endl;\
            return 1;\
        }\
        if(options.count_options("help")>0){\
            string help = options.get_value<string>("help","");\
            if(help=="traj"){\
                cout << engine.help() << endl;\
            } else if(help=="task"){\
                cout << task.help() << endl;\
            } else {\
                cout << task.help() << endl;\
                cout << engine.help() << endl;\
            }\
            return 1;\
        }\
        engine.run();\
    } catch(const Pteros_error& e) {\
        cout << e.what() << endl;\
    }\
}


#endif //STANDALONE_PLUGINS

#endif //COMPILED_PLUGIN_H
