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

#pragma once

#include "pteros/analysis/task_plugin.h"

#ifndef STANDALONE_PLUGINS

// Make a Python extension module from this plugin

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>

namespace py = pybind11;

#define _EVAL(arg) #arg
#define CREATE_COMPILED_PLUGIN(_name) \
PYBIND11_MODULE(_name, m) {\
    py::class_<_name,TaskPlugin,std::shared_ptr<_name>>(m, _EVAL(_name))\
        .def(py::init<const Options&>())\
        .def("help",&_name::help)\
    ;\
}

#else //STANDALONE_PLUGINS

// Make a stand-alone executable from this plugin

#include "pteros/analysis/trajectory_reader.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

using namespace pteros;
using namespace std;

// Skeleton of the driver program

#define CREATE_COMPILED_PLUGIN(_name)\
int main(int argc, char** argv){\
    try {\
        Options options;\
        parse_command_line(argc,argv,options);\
        Trajectory_reader engine(options);\
        auto task = new _name(options);\
        engine.add_task(task);\
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




