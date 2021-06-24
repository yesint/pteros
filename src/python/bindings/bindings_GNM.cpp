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


#include "bindings_util.h"
#include "pteros/extras/gnm.h"

namespace py = pybind11;
using namespace pteros;
using namespace pybind11::literals;

void make_bindings_GNM(py::module& m){

    py::class_<GNM>(m,"GNM")
        .def(py::init<const Selection&,float>())
        .def(py::init<const Selection&>())
        .def("get_eigenvector",&GNM::get_eigenvector)
        .def("get_B_factor",&GNM::get_B_factor)
        .def("get_subset_c_matrix",&GNM::get_subset_c_matrix)
        .def("get_c_matrix",&GNM::get_c_matrix)
        .def("write_eigenvectors",&GNM::write_eigenvectors)
        .def("compute_c_matrix",&GNM::compute_c_matrix)
        .def("write_c_matrix",&GNM::write_c_matrix)
    ;

}




