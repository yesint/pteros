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

#include "bindings_options.h"
#include "pteros/python/bindings_util.h"
#include "pteros/analysis/options.h"

using namespace pteros;
using namespace Eigen;
using namespace boost::python;

const Option& Options_call1(Options* o, std::string key){
    return (*o)(key);
}

const Option& Options_call2(Options* o, std::string key, std::string default_val){
    return (*o)(key,default_val);
}

// Version for parsing python argv
Options parse_command_line1(boost::python::list py_argv){
    int argc = len(py_argv);
    char* argv[argc];
    string tmp;
    // Allocate memory for c-strings and copy
    for(int i=0;i<argc;++i){
        argv[i] = (char*)malloc( sizeof(char) * len(py_argv[i]) );
        tmp = extract<string>(py_argv[i]);
        strcpy(argv[i], tmp.c_str());
    }
    // Call
    Options opt;
    parse_command_line(argc,argv,opt);
    // Free memory of c-strings
    for(int i=0;i<argc;++i) free(argv[i]);

    return opt;
}

boost::python::list Option_as_ints(Option* o){
    boost::python::list res;
    auto vec = o->as_ints();
    for(auto& v: vec) res.append(v);
    return res;
}

boost::python::list Option_as_floats(Option* o){
    boost::python::list res;
    auto vec = o->as_floats();
    for(auto& v: vec) res.append(v);
    return res;
}

boost::python::list Option_as_bools(Option* o){
    boost::python::list res;
    auto vec = o->as_bools();
    for(bool v: vec) res.append(v);
    return res;
}

boost::python::list Option_as_strings(Option* o){
    boost::python::list res;
    auto vec = o->as_strings();
    for(auto& v: vec) res.append(v);
    return res;
}


void make_bindings_Options(){

    import_array();

    def("parse_command_line", &parse_command_line1);

    class_<Option>("Option", init<>())
        .def("as_int",&Option::as_int)
        .def("as_float",&Option::as_float)
        .def("as_bool",&Option::as_bool)
        .def("as_string",&Option::as_string)

        .def("as_ints",&Option_as_ints)
        .def("as_floats",&Option_as_floats)
        .def("as_bools",&Option_as_bools)
        .def("as_strings",&Option_as_strings)
    ;

    class_<Options>("Options", init<>())
        .def("__call__",&Options_call1, return_value_policy<reference_existing_object>())
        .def("__call__",&Options_call2, return_value_policy<reference_existing_object>())
        .def("debug",&Options::debug)
    ;
}
