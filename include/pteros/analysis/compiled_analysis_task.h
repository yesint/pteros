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

#ifndef COMPILED_ANALYSIS_TASK_H
#define COMPILED_ANALYSIS_TASK_H

#include "pteros/analysis/trajectory_reader.h"
#include "pteros/analysis/task_plugin.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

using namespace pteros;
using namespace std;

// Skeleton of the driver program

#define COMPILED_ANALYSIS_TASK(_name, _collector)\
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


#endif
