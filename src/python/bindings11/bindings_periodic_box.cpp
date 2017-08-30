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

#include "pteros/core/periodic_box.h"
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace Eigen;
using namespace pybind11::literals;

void make_bindings_Periodic_box(py::module& m){

    py::class_<Periodic_box>(m, "Periodic_box")
        .def(py::init<Vector3f_const_ref,Vector3f_const_ref>())
        .def("get_vector",&Periodic_box::get_vector)
        .def("get_matrix",[](Periodic_box* b){Matrix3f m = b->get_matrix().transpose(); return m;})
        .def("set_matrix",[](Periodic_box* b, Matrix3f_const_ref m){ b->set_matrix(m.transpose()); })
        .def("scale_vectors",&Periodic_box::scale_vectors)
        .def("get_inv_matrix",[](Periodic_box* b){Matrix3f m = b->get_inv_matrix().transpose(); return m;})
        .def("lab_to_box",&Periodic_box::lab_to_box)
        .def("lab_to_box_transform",[](Periodic_box* b){Matrix3f m = b->lab_to_box_transform().transpose(); return m;})
        .def("box_to_lab",&Periodic_box::box_to_lab)
        .def("box_to_lab_transform",[](Periodic_box* b){Matrix3f m = b->box_to_lab_transform().transpose(); return m;})
        .def("extent",&Periodic_box::extent)
        .def("extents",&Periodic_box::extents)
        .def("is_triclinic",&Periodic_box::is_triclinic)
        .def("is_periodic",&Periodic_box::is_periodic)
        .def("distance",&Periodic_box::distance, "point1"_a, "point2"_a, "dims"_a=Eigen::Vector3i::Ones())
        .def("distance_squared",&Periodic_box::distance_squared, "point1"_a, "point2"_a, "dims"_a=Eigen::Vector3i::Ones())
        .def("wrap_point",&Periodic_box::wrap_point, "point"_a, "dims"_a=Eigen::Vector3i::Ones())
        .def("in_box",&Periodic_box::in_box, "point"_a, "origin"_a=Eigen::Vector3i::Zero())
        .def("get_closest_image",&Periodic_box::get_closest_image, "point"_a, "target"_a, "dims"_a=Eigen::Vector3i::Ones())
        .def("shortest_vector",&Periodic_box::shortest_vector, "point1"_a, "point2"_a, "dims"_a=Eigen::Vector3i::Ones())
        .def("volume",&Periodic_box::volume)
        .def("read_pdb_box",&Periodic_box::read_pdb_box)
        .def("write_pdb_box",&Periodic_box::write_pdb_box)
        .def("from_vectors_angles",&Periodic_box::from_vectors_angles)
        .def("to_vectors_angles",[](Periodic_box* b){
            Vector3f vec;
            Vector3f ang;
            b->to_vectors_angles(vec,ang);
            return py::make_tuple(vec,ang);
        })
    ;
}
