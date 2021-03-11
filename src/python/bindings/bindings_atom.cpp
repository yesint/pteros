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


#include "pteros/core/atom.h"
#include "bindings_util.h"

namespace py = pybind11;
using namespace pteros;

void make_bindings_Atom(py::module& m){

    py::class_<Atom>(m, "Atom")
        .def(py::init<>())
        .def_readwrite("resid",     &Atom::resid)
        .def_readwrite("name",      &Atom::name)
        .def_readwrite("chain",     &Atom::chain)
        .def_readwrite("resname",   &Atom::resname)
        .def_readwrite("tag",       &Atom::tag)
        .def_readwrite("occupancy", &Atom::occupancy)
        .def_readwrite("beta",      &Atom::beta)
        .def_readwrite("resindex",  &Atom::resindex)
        .def_readwrite("mass",      &Atom::mass)
        .def_readwrite("charge",    &Atom::charge)
        .def_readwrite("type",      &Atom::type)
        .def_readwrite("type_name", &Atom::type_name)
        .def_readwrite("atomic_number", &Atom::atomic_number)
    ;
}




