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


#include "pteros/analysis/options.h"
#include "bindings_util.h"


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




