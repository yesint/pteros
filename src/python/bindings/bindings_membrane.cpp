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

#include "bindings_util.h"
#include "pteros/extras/membrane.h"

namespace py = pybind11;
using namespace pteros;
using namespace pybind11::literals;

void make_bindings_Membrane(py::module& m){

    py::class_<Lipid_descr>(m,"Lipid_descr")
        .def(py::init<const std::string&,const std::string&,const std::string&,const std::string&,const std::string&>())
        .def_readonly("name",&Lipid_descr::name)
        .def_readonly("whole_sel_str",&Lipid_descr::whole_sel_str)
        .def_readonly("head_sel_str",&Lipid_descr::head_sel_str)
        .def_readonly("tail_sel_str",&Lipid_descr::tail_sel_str)
        .def_readonly("head_mid_str",&Lipid_descr::mid_sel_str)
    ;

    py::class_<Splay_pair>(m,"Splay_pair")
        .def_readonly("lip1",&Splay_pair::lip1)
        .def_readonly("lip2",&Splay_pair::lip2)
        .def_readonly("splay",&Splay_pair::splay)
    ;

    py::class_<Lipid>(m,"Lipid")
            .def(py::init<const Selection&,const Lipid_descr&>())

            .def("get_mid_xyz",&Lipid::get_mid_xyz)
            .def("get_head_xyz",&Lipid::get_head_xyz)
            .def("get_tail_xyz",&Lipid::get_tail_xyz)

            .def_readwrite("mid_sel",&Lipid::mid_sel)
            .def_readwrite("head_sel",&Lipid::head_sel)
            .def_readwrite("tail_sel",&Lipid::tail_sel)
            .def_readwrite("whole_sel",&Lipid::whole_sel)

            .def_readonly("name",&Lipid::name)
            .def_readonly("normal",&Lipid::normal)
            .def_readonly("smoothed_mid_xyz",&Lipid::smoothed_mid_xyz)
            .def_readonly("tilt",&Lipid::tilt)            
            .def_readonly("quad_fit_rms",&Lipid::quad_fit_rms)
            .def_readonly("area",&Lipid::area)
            .def_readonly("leaflet",&Lipid::leaflet)
            .def_readonly("gaussian_curvature",&Lipid::gaussian_curvature)
            .def_readonly("mean_curvature",&Lipid::mean_curvature)
            .def_readonly("coord_number",&Lipid::coord_number)
    ;

    py::class_<Membrane>(m,"Membrane")
        .def(py::init<System*,const std::vector<Lipid_descr>&>())
        .def("compute_properties",&Membrane::compute_properties,"d"_a,"external_normal"_a=Eigen::Vector3f::Zero())
        .def("write_vmd_arrows",&Membrane::write_vmd_arrows)        
        .def("write_smoothed",&Membrane::write_smoothed)
        .def("num_lipids",&Membrane::num_lipids)
        .def("get_lipid",&Membrane::get_lipid,py::return_value_policy::reference_internal)
        .def("num_leaflets",&Membrane::num_leaflets)
        .def("get_leaflet",&Membrane::get_leaflet)
        .def_readonly("lipids",&Membrane::lipids)
        .def_readonly("splay",&Membrane::splay)
        .def_readonly("neighbors",&Membrane::neighbors)
        .def_readonly("leaflets",&Membrane::leaflets)
    ;
}
