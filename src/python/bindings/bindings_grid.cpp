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


#include "pteros/core/grid.h"
#include "bindings_util.h"

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace pybind11::literals;

void make_bindings_Grid(py::module& m){

    py::class_<GridCell>(m, "GridCell")
        .def("add_point",&GridCell::add_point)
        .def("clear",&GridCell::clear)
        .def("get_index",&GridCell::get_index)
        .def("get_coord",&GridCell::get_coord)
        .def("size",&GridCell::size)
        .def("get_average_coord",&GridCell::get_average_coord)
    ;

    py::class_<Grid>(m, "Grid")
        .def(py::init<>())
        .def(py::init<int,int,int>())
        .def("cell",py::overload_cast<int,int,int>(&Grid::cell))
        .def("cell",py::overload_cast<Vector3i_const_ref>(&Grid::cell))
        .def("resize",&Grid::resize)
        .def("clear",&Grid::clear)
        .def("populate",py::overload_cast<const Selection&,bool>(&Grid::populate),"sel"_a,"abs_index"_a=false)
        .def("populate",py::overload_cast<const Selection&,Vector3f_const_ref,Vector3f_const_ref,bool>(&Grid::populate),"sel"_a,"min"_a,"max"_a,"abs_index"_a=false)
        .def("populate_periodic",py::overload_cast<const Selection&,Vector3i_const_ref,bool>(&Grid::populate_periodic),"sel"_a,"pbc_dims"_a=fullPBC,"abs_index"_a=false)
        .def("populate_periodic",py::overload_cast<const Selection&,const PeriodicBox&,Vector3i_const_ref,bool>(&Grid::populate_periodic),"sel"_a,"box"_a,"pbc_dims"_a=fullPBC,"abs_index"_a=false)
    ;
}




