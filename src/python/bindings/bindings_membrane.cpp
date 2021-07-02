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
#include "pteros/extras/membrane.h"

namespace py = pybind11;
using namespace pteros;
using namespace pybind11::literals;

void make_bindings_Membrane(py::module& m){

    py::class_<LipidSpecies>(m,"LipidSpecies")
        .def(py::init<const std::string&,
             const std::string&,
             const std::string&,
             const std::string&,
             const std::string&,
             const std::vector<std::string>&>())
        .def_readonly("name",&LipidSpecies::name)
        .def_readonly("whole_str",&LipidSpecies::whole_str)
        .def_readonly("head_marker_str",&LipidSpecies::head_marker_str)
        .def_readonly("tail_marker_str",&LipidSpecies::tail_marker_str)
        .def_readonly("mid_marker_str",&LipidSpecies::mid_marker_str)
        .def_readonly("tail_carbons_str",&LipidSpecies::tail_carbons_str)
    ;


    py::class_<LipidMolecule>(m,"LipidMolecule")
        //.def(py::init<const Selection&,const LipidSpecies&,int,LipidMembrane*>())
        .def_readonly("name",&LipidMolecule::name)
        .def_readwrite("mid_marker_sel",&LipidMolecule::mid_marker_sel)
        .def_readwrite("head_marker_sel",&LipidMolecule::head_marker_sel)
        .def_readwrite("tail_marker_sel",&LipidMolecule::tail_marker_sel)
        .def_readwrite("whole_sel",&LipidMolecule::whole_sel)

        .def("add_to_group",&LipidMolecule::add_to_group)
    ;

    py::class_<LipidMembrane>(m,"LipidMembrane")
        .def(py::init<System*,const std::vector<LipidSpecies>&,int>())

        .def("compute_properties",&LipidMembrane::compute_properties,
             "d"_a=2.0,
             "use_external_normal"_a=false,
             "external_pivot"_a=Eigen::Vector3f::Zero(),
             "external_dist_dim"_a=Eigen::Vector3i::Ones()
            )
        .def("reset_groups",&LipidMembrane::reset_groups)
        .def("compute_averages",&LipidMembrane::compute_averages)
        .def("write_averages",&LipidMembrane::write_averages,"path"_a="")

        .def_readonly("lipids",&LipidMembrane::lipids)
    ;
}




