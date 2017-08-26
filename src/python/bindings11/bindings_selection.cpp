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
#include <pybind11/operators.h>

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace pybind11::literals;

void make_bindings_Selection(py::module& m){

    using RowMatrixXf = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    py::class_<Selection>(m, "Selection")
        // Constructors
        .def(py::init<>())
        .def(py::init<const System&>())
        .def(py::init<const System&,std::string,int>(),"sys"_a,"str"_a,"fr"_a=0)
        .def(py::init<const System&,int,int>())
        .def(py::init<const System&,const std::vector<int>&>())
        .def(py::init<const System&,const std::function<void(const System&,int,std::vector<int>&)>&,int>(),"sys"_a,"callback"_a,"fr"_a=0)
        .def(py::init<const Selection&>())

        // Oparators
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self | py::self)
        .def(py::self & py::self)
        .def(py::self - py::self)
        .def(~py::self)

        // Modification
        .def("append", py::overload_cast<const Selection&>(&Selection::append))
        .def("append", py::overload_cast<int>(&Selection::append))
        .def("remove", py::overload_cast<const Selection&>(&Selection::remove))
        .def("remove", py::overload_cast<int>(&Selection::remove))
        .def("invert",&Selection::invert)
        .def("set_system",&Selection::set_system)
        .def("modify", py::overload_cast<string,int>(&Selection::modify),"str"_a,"fr"_a=0)
        .def("modify", py::overload_cast<int,int>(&Selection::modify))
        .def("modify", py::overload_cast<const std::vector<int>&>(&Selection::modify))
        .def("modify", py::overload_cast<const std::function<void(const System&,int,std::vector<int>&)>&,int>(&Selection::modify),"callback"_a,"fr"_a=0)
        .def("modify", py::overload_cast<const System&,string,int>(&Selection::modify),"sys"_a,"str"_a,"fr"_a=0)
        .def("modify", py::overload_cast<const System&,int,int>(&Selection::modify))
        .def("modify", py::overload_cast<const System&,const std::vector<int>&>(&Selection::modify))
        .def("modify", py::overload_cast<const System&,const std::function<void(const System&,int,std::vector<int>&)>&,int>(&Selection::modify),"sys"_a,"callback"_a,"fr"_a=0)
        .def("apply",&Selection::apply)
        .def("update",&Selection::update)
        .def("clear",&Selection::clear)

        // Subselection
        .def("__call__", py::overload_cast<string>(&Selection::operator()))
        .def("__call__", py::overload_cast<int,int>(&Selection::operator()))
        .def("__call__", py::overload_cast<const std::vector<int>&>(&Selection::operator()))
        .def("select", py::overload_cast<string>(&Selection::operator()))
        .def("select", py::overload_cast<int,int>(&Selection::operator()))
        .def("select", py::overload_cast<const std::vector<int>&>(&Selection::operator()))

         .def("size",&Selection::size)
        .def("__len__", &Selection::size)
        .def("get_index", &Selection::get_index)
        .def("translate", &Selection::translate)
        .def("get_xyz", [](Selection* sel){  RowMatrixXf m = sel->get_xyz().transpose(); return m; })
        .def("set_xyz", [](Selection* sel, MatrixXf_const_ref m){ sel->set_xyz(m.transpose()); })
    ;
}
