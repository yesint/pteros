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

#include "pteros/analysis/options.h"
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>


namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace pybind11::literals;

void make_bindings_Options(py::module& m){

    py::class_<Option>(m, "Option")
        .def("as_string",&Option::as_string)
        .def("as_int",&Option::as_int)
        .def("as_float",&Option::as_float)
        .def("as_bool",&Option::as_bool)

        .def("as_strings",&Option::as_strings)
        .def("as_ints",&Option::as_ints)
        .def("as_floats",&Option::as_floats)
        // vector of bools fails to convert automatically in pybind11...
        .def("as_bools",[](Option* o){
            py::list res;
            auto vec = o->as_bools();
            for(bool v: vec) res.append(v);
            return res;
        })
    ;

    py::class_<Options>(m, "Options")
        .def("__call__", py::overload_cast<string>(&Options::operator(), py::const_), "key"_a)
        .def("__call__", py::overload_cast<string,string>(&Options::operator()), "key"_a, "default_value"_a)
        .def("has",&Options::has)
        .def("get_name",&Options::get_name)
    ;

    // parser with tasks
    m.def("parse_command_line", [](const py::list& py_argv, const std::string& task_tag){
        // Make c-style argc/argv from python argv
        int argc = py::len(py_argv);
        char* argv[argc];
        string tmp;
        // Allocate memory for c-strings and copy
        for(int i=0;i<argc;++i){
            argv[i] = (char*)malloc( sizeof(char) * py::len(py_argv[i]) );
            tmp = py_argv[i].cast<string>();
            strcpy(argv[i], tmp.c_str());
        }
        // Call
        Options opt;
        std::vector<Options> tasks;
        parse_command_line(argc,argv,opt,task_tag,tasks);
        // Free memory of c-strings
        for(int i=0;i<argc;++i) free(argv[i]);
        // return
        return py::make_tuple(opt,tasks);
    }, "argv"_a, "task_tag"_a);

    // parser without tasks
    m.def("parse_command_line", [](const py::list& py_argv){
        // Make c-style argc/argv from python argv
        int argc = py::len(py_argv);
        char* argv[argc];
        string tmp;
        // Allocate memory for c-strings and copy
        for(int i=0;i<argc;++i){
            argv[i] = (char*)malloc( sizeof(char) * py::len(py_argv[i]) );
            tmp = py_argv[i].cast<string>();
            strcpy(argv[i], tmp.c_str());
        }
        // Call
        Options opt;
        parse_command_line(argc,argv,opt);
        // Free memory of c-strings
        for(int i=0;i<argc;++i) free(argv[i]);
        // return
        return opt;
    }, "argv"_a);

}
