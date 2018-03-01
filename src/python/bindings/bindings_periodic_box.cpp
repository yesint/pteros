/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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


#include "pteros/core/periodic_box.h"
#include "bindings_util.h"

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace Eigen;
using namespace pybind11::literals;

void make_bindings_Periodic_box(py::module& m){

    py::class_<Periodic_box>(m, "Periodic_box")
        .def(py::init<Vector3f_const_ref,Vector3f_const_ref>())
        .def(py::init<const Periodic_box&>())
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
        .def("distance",&Periodic_box::distance, "point1"_a, "point2"_a, "pbc"_a=fullPBC)
        .def("distance_squared",&Periodic_box::distance_squared, "point1"_a, "point2"_a, "pbc"_a=fullPBC)
        .def("wrap_point",&Periodic_box::wrap_point, "point"_a, "pbc"_a=fullPBC, "origin"_a=Eigen::Vector3f::Zero())
        .def("in_box",&Periodic_box::in_box, "point"_a, "origin"_a=noPBC)
        .def("closest_image",&Periodic_box::closest_image, "point"_a, "target"_a, "pbc"_a=fullPBC)
        .def("shortest_vector",&Periodic_box::shortest_vector, "point1"_a, "point2"_a, "pbc"_a=fullPBC)
        .def("volume",&Periodic_box::volume)
        .def("from_pdb_box",&Periodic_box::from_pdb_box)
        .def("to_pdb_box",&Periodic_box::to_pdb_box)
        .def("from_vectors_angles",&Periodic_box::from_vectors_angles)
        .def("to_vectors_angles",[](Periodic_box* b){
            Vector3f vec;
            Vector3f ang;
            b->to_vectors_angles(vec,ang);
            return py::make_tuple(vec,ang);
        })
    ;
}

