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


#include "pteros/core/periodic_box.h"
#include "bindings_util.h"

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace Eigen;
using namespace pybind11::literals;

void make_bindings_Periodic_box(py::module& m){

    py::class_<PeriodicBox>(m, "PeriodicBox")
        .def(py::init<>())
        .def(py::init<Vector3f_const_ref,Vector3f_const_ref>())
        .def(py::init<const PeriodicBox&>())
        .def("get_element",&PeriodicBox::get_element)
        .def("set_element",&PeriodicBox::set_element)
        .def("get_vector",&PeriodicBox::get_vector)
        .def("set_vector",&PeriodicBox::set_vector)
        .def("get_matrix",[](PeriodicBox* b){Matrix3f m = b->get_matrix().transpose(); return m;})
        .def("set_matrix",[](PeriodicBox* b, Matrix3f_const_ref m){ b->set_matrix(m.transpose()); })
        .def("scale_vectors",&PeriodicBox::scale_vectors)
        .def("get_inv_matrix",[](PeriodicBox* b){Matrix3f m = b->get_inv_matrix().transpose(); return m;})
        .def("lab_to_box",&PeriodicBox::lab_to_box)
        .def("lab_to_box_transform",[](PeriodicBox* b){Matrix3f m = b->lab_to_box_transform().transpose(); return m;})
        .def("box_to_lab",&PeriodicBox::box_to_lab)
        .def("box_to_lab_transform",[](PeriodicBox* b){Matrix3f m = b->box_to_lab_transform().transpose(); return m;})
        .def("extent",&PeriodicBox::extent)
        .def("extents",&PeriodicBox::extents)
        .def("is_triclinic",&PeriodicBox::is_triclinic)
        .def("is_periodic",&PeriodicBox::is_periodic)
        .def("distance",&PeriodicBox::distance, "point1"_a, "point2"_a, "pbc"_a=fullPBC)
        .def("distance_squared",&PeriodicBox::distance_squared, "point1"_a, "point2"_a, "pbc"_a=fullPBC)
        .def("wrap_point",&PeriodicBox::wrap_point, "point"_a, "pbc"_a=fullPBC, "origin"_a=Eigen::Vector3f::Zero())
        .def("in_box",&PeriodicBox::in_box, "point"_a, "origin"_a=noPBC)
        .def("closest_image",&PeriodicBox::closest_image, "point"_a, "target"_a, "pbc"_a=fullPBC)
        .def("shortest_vector",&PeriodicBox::shortest_vector, "point1"_a, "point2"_a, "pbc"_a=fullPBC)
        .def("volume",&PeriodicBox::volume)
        .def("from_pdb_box",&PeriodicBox::from_pdb_box)
        .def("to_pdb_box",&PeriodicBox::to_pdb_box)
        .def("from_vectors_angles",&PeriodicBox::from_vectors_angles)
        .def("to_vectors_angles",[](PeriodicBox* b){
            Vector3f vec;
            Vector3f ang;
            b->to_vectors_angles(vec,ang);
            return py::make_tuple(vec,ang);
        })
    ;
}




