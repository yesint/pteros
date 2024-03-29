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


#include "bindings_util.h"
#include "pteros/extras/membrane/lipid_membrane.h"

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
        .def_readonly("whole_sel_str",&LipidSpecies::whole_sel_str)
        .def_readonly("head_sel_names",&LipidSpecies::head_subsel_names)
        .def_readonly("tail_sel_names",&LipidSpecies::tail_subsel_names)
        .def_readonly("surf_marker_names",&LipidSpecies::surf_subsel_names)

        //.def_readonly("tail_carbons_str",&LipidSpecies::tails_descr)
    ;


    py::class_<LipidMolecule>(m,"LipidMolecule")
        //.def(py::init<const Selection&,const LipidSpecies&,int,LipidMembrane*>())
        //.def_readonly("name",&LipidMolecule::name)
        .def_readwrite("mid_marker_sel",&LipidMolecule::surf_marker_sel)
        .def_readwrite("head_marker_sel",&LipidMolecule::head_marker_sel)
        .def_readwrite("tail_marker_sel",&LipidMolecule::tail_marker_sel)
        .def_readwrite("whole_sel",&LipidMolecule::whole_sel)

        .def("add_to_group",&LipidMolecule::add_to_group)

        .def_readonly("normal",&LipidMolecule::normal)
        .def_readonly("tilt",&LipidMolecule::tilt)
        .def_readonly("mean_curvature",&LipidMolecule::mean_curvature)
        .def_readonly("gaussian_curvature",&LipidMolecule::gaussian_curvature)

        .def_readonly("coord_number",&LipidMolecule::coord_number)

        .def_readonly("neib",&LipidMolecule::neib)
        .def_readonly("smoothed_mid_xyz",&LipidMolecule::smoothed_surf_marker)
        .def_readonly("quad_fit_rms",&LipidMolecule::quad_fit_rms)

    ;

    py::enum_<OrderType>(m,"OrderType")
        .value("SZ",OrderType::SZ)
        .value("SCD",OrderType::SCD)
        .value("SCD_CORR",OrderType::SCD_CORR)
    ;

    py::class_<InterpolatedPoint>(m,"InterpolatedPoint")
        .def_readonly("normal",&InterpolatedPoint::normal)
        .def_readonly("neib_lipids",&InterpolatedPoint::neib_lipids)
        .def_readonly("weights",&InterpolatedPoint::weights)
        .def_readonly("mean_curvature",&InterpolatedPoint::mean_curvature)
        .def_readonly("mean_depth",&InterpolatedPoint::mean_depth)
    ;

    py::class_<LipidMembraneOptions>(m,"LipidMembraneOptions")
        .def_readwrite("smoothing_tol",&LipidMembraneOptions::smoothing_tol)
        .def_readwrite("smoothing_maxiter",&LipidMembraneOptions::smoothing_maxiter)
        .def_readwrite("per_carbon_normals",&LipidMembraneOptions::per_carbon_normals)
        .def_readwrite("inclusion_h_cutoff",&LipidMembraneOptions::inclusion_h_cutoff)
        .def("__enter__",[](LipidMembraneOptions* self)->LipidMembraneOptions* {return self;})
        .def("__exit__",[](LipidMembraneOptions* self,
             const py::object &type,
             const py::object &value,
             const py::object &traceback){})
    ;

    py::class_<LipidMembrane>(m,"LipidMembrane")
        .def(py::init<const Selection&,int,const std::vector<LipidSpecies>&,const Selection&>(),
             "input_sel"_a,"ngroups"_a,"sp_list"_a,"incl"_a=Selection())

        .def_static("get_domains",&LipidMembrane::get_domains,"sys"_a,"sp_list"_a,"d"_a=0.4)
        .def("compute_properties",&LipidMembrane::compute_properties,
             "d"_a=2.0, "incl_d"_a=0.5, "order_type"_a=OrderType::SCD_CORR)
        .def("reset_groups",&LipidMembrane::reset_groups)
        .def("compute_averages",&LipidMembrane::compute_averages)
        .def("get_average_curvatures",&LipidMembrane::get_average_curvatures)
        .def("compute_triangulation",&LipidMembrane::compute_triangulation)
        .def("write_averages",&LipidMembrane::write_averages,"path"_a=".")
        .def("write_vmd_visualization",&LipidMembrane::write_vmd_visualization,"path"_a=".")
        .def("get_interpolation",[](LipidMembrane* lm, const Selection& points, float d){
            std::vector<InterpolatedPoint> res;
            Eigen::MatrixXf m;
            lm->get_interpolation(points,res,m,d);
            return py::make_tuple(res,m);
        }, "points"_a, "d"_a=3.0)

        .def_readonly("lipids",&LipidMembrane::lipids)
        .def_readwrite("options",&LipidMembrane::options)
    ;
}




