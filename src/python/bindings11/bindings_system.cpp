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
using namespace Eigen;
using namespace pybind11::literals;

void make_bindings_System(py::module& m){

    py::class_<System>(m, "System")
        .def(py::init<const std::string &>())

        // Size
        .def("num_atoms", &System::num_atoms)
        .def("num_frames", &System::num_frames)

        // Append
        .def("append", py::overload_cast<const System&>(&System::append))
        .def("append", py::overload_cast<const Selection&>(&System::append))
        .def("append", py::overload_cast<const Atom&, Vector3f_const_ref>(&System::append))

        // Reaaranging
        .def("rearrange", py::overload_cast<const vector<string>&>(&System::rearrange))
        .def("rearrange", py::overload_cast<const vector<Selection>&>(&System::rearrange))
        .def("keep", py::overload_cast<const string&>(&System::keep))
        .def("keep", py::overload_cast<const Selection&>(&System::keep))
        .def("remove", py::overload_cast<const string&>(&System::remove))
        .def("remove", py::overload_cast<Selection&>(&System::remove))
        .def("distribute", [](System* s, const Selection sel, Vector3i_const_ref ncopies, Matrix3f_const_ref shift){
            Matrix3f m = shift.transpose();
            s->distribute(sel,ncopies,m);
        })

        // Loading
        .def("load", py::overload_cast<string,int,int,int,std::function<bool(System*,int)>>(&System::load),
             "fname"_a, "b"_a=0, "e"_a=-1, "skip"_a=0, "on_frame"_a=nullptr)

        // Selecting
        .def("__call__", py::overload_cast<>(&System::operator()))
        .def("__call__", py::overload_cast<string,int>(&System::operator()),"str"_a,"fr"_a=0)
        .def("__call__", py::overload_cast<int,int>(&System::operator()))
        .def("__call__", py::overload_cast<const std::vector<int>&>(&System::operator()))
        .def("__call__", py::overload_cast<const std::function<void(const System&,int,std::vector<int>&)>&,int>(&System::operator()),"callback"_a,"fr"_a=0)
        .def("select_all", py::overload_cast<>(&System::operator()))
        .def("select", py::overload_cast<string,int>(&System::operator()),"str"_a,"fr"_a=0)
        .def("select", py::overload_cast<int,int>(&System::operator()))
        .def("select", py::overload_cast<const std::vector<int>&>(&System::operator()))
        .def("select", py::overload_cast<const std::function<void(const System&,int,std::vector<int>&)>&,int>(&System::operator()),"callback"_a,"fr"_a=0)

        // Input filtering
        .def("set_filter", py::overload_cast<string>(&System::set_filter))
        .def("set_filter", py::overload_cast<int,int>(&System::set_filter))
        .def("set_filter", py::overload_cast<const std::vector<int>&>(&System::set_filter))

        // Frame operations
        .def("frame_dup", &System::frame_dup)
        .def("frame_append", &System::frame_append)
        .def("frame_copy", &System::frame_copy)
        .def("frame_delete", &System::frame_delete, "b"_a=0, "e"_a=-1)
        .def("frame_swap", &System::frame_swap)

        // Accessors
        .def("getBox", py::overload_cast<int>(&System::Box, py::const_))
        .def("setBox", [](System* s,int fr,const Periodic_box& b){ s->Box(fr)=b; })

        .def("getTime", py::overload_cast<int>(&System::Time, py::const_))
        .def("setTime", [](System* s,int fr,int t){ s->Time(fr)=t; })

        .def("getFrame_data", py::overload_cast<int>(&System::Frame_data, py::const_))
        .def("setFrame_data", [](System* s, int i, const Frame& fr){ s->Frame_data(i)=fr; })

        .def("getXYZ", py::overload_cast<int,int>(&System::XYZ, py::const_))
        .def("setXYZ", [](System* s,Vector3f_const_ref v,int i,int fr){ s->XYZ(i,fr)=v; })

        .def("getAtom_data", py::overload_cast<int>(&System::Atom_data, py::const_))
        .def("setAtom_data", [](System* s, int i, const Atom& a){ s->Atom_data(i)=a; })

        // dssp
        .def("dssp", py::overload_cast<string,int>(&System::dssp, py::const_))
        .def("dssp", py::overload_cast<int>(&System::dssp, py::const_))

        // operations with atoms
        .def("atoms_dup", &System::atoms_dup)
        .def("atoms_add", &System::atoms_add)
        .def("atoms_delete", &System::atoms_delete)
        .def("atom_move", &System::atom_move)
        .def("atom_clone", &System::atom_clone)

        // wrap
        .def("wrap", &System::wrap, "fr"_a, "dims"_a=Eigen::Vector3i::Ones())

        // measuring
        .def("distance", &System::distance, "i"_a, "j"_a, "fr"_a, "periodic"_a=true, "dims"_a=Eigen::Vector3i::Ones())
        .def("angle", &System::angle, "i"_a, "j"_a, "k"_a, "fr"_a, "periodic"_a=true, "dims"_a=Eigen::Vector3i::Ones())
        .def("dihedral", &System::dihedral, "i"_a, "j"_a, "k"_a, "l"_a, "fr"_a, "periodic"_a=true, "dims"_a=Eigen::Vector3i::Ones())

        // Energy
        .def("non_bond_energy", py::overload_cast<int,int,int,bool>(&System::non_bond_energy,py::const_), "at1"_a, "at2"_a, "fr"_a, "periodic"_a=true)
        .def("non_bond_energy", py::overload_cast<const std::vector<Eigen::Vector2i>&,int,bool>(&System::non_bond_energy,py::const_), "nlist"_a, "fr"_a, "periodic"_a=true)

        // Unit
        .def("clear", &System::clear)
        .def("force_field_ready", &System::force_field_ready)
        .def("assign_resindex", &System::assign_resindex, "start"_a=0)
        .def("sort_by_resindex", &System::sort_by_resindex)

    ;
}
