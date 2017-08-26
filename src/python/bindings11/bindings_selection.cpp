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

#include "pteros/core/selection.h"
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace pybind11::literals;

void make_bindings_Selection(py::module& m){

    using RowMatrixXf = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    py::class_<Selection>(m, "Selection")
        .def("size",&Selection::size)
        .def("__len__", &Selection::size)
        .def("get_index", &Selection::get_index)
        .def("translate", &Selection::translate)
        .def("get_xyz", [](Selection* sel){  RowMatrixXf m = sel->get_xyz().transpose(); return m; })
        .def("set_xyz", [](Selection* sel, MatrixXf_const_ref m){ sel->set_xyz(m.transpose()); })
    ;
}
