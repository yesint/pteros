/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2023, Semen Yesylevskyy
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


#include "pteros/core/system.h"
#include "bindings_util.h"

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace pybind11::literals;

void make_bindings_Frame(py::module& m){

    py::class_<Frame>(m, "Frame")
        .def_readwrite("time", &Frame::time)
        .def_readwrite("box", &Frame::box)
        .def_readwrite("coord", &Frame::coord)
        .def("has_vel",&Frame::has_vel)
        .def("has_force",&Frame::has_force)
    ;
}




