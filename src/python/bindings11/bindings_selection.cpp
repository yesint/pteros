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
#include "pteros/core/pteros_error.h"

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace Eigen;
using namespace pybind11::literals;

#define DEF_PROPERTY(_name,_dtype,_func) \
    .def_property(#_name, [](Atom_proxy* obj){return obj->_func();}, [](Atom_proxy* obj,const _dtype& val){obj->_func()=val;})

void make_bindings_Selection(py::module& m){

    using RowMatrixXf = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    py::class_<Atom_proxy>(m, "Atom_proxy")
        .def_property_readonly("index", [](Atom_proxy* obj){return obj->Index();})
        DEF_PROPERTY(resid,int,Resid)
        DEF_PROPERTY(resindex,int,Resindex)
        DEF_PROPERTY(resname,string,Resname)
        DEF_PROPERTY(name,string,Name)
        DEF_PROPERTY(chain,char,Chain)
        DEF_PROPERTY(tag,string,Tag)
        DEF_PROPERTY(occupancy,float,Occupancy)
        DEF_PROPERTY(beta,float,Beta)
        DEF_PROPERTY(mass,float,Mass)
        DEF_PROPERTY(charge,float,Charge)
        DEF_PROPERTY(type,int,Type)
        DEF_PROPERTY(type_name,string,Type_name)
        DEF_PROPERTY(atom,Atom,Atom_data)
        DEF_PROPERTY(x,float,X)
        DEF_PROPERTY(y,float,Y)
        DEF_PROPERTY(z,float,Z)
        .def_property("xyz", [](Atom_proxy* obj){return obj->XYZ();}, [](Atom_proxy* obj,Vector3f_const_ref val){obj->XYZ()=val;})
        // For other frame
        .def("getX", [](Atom_proxy* ap, int fr){ return ap->X(fr); })
        .def("getY", [](Atom_proxy* ap, int fr){ return ap->Y(fr); })
        .def("getZ", [](Atom_proxy* ap, int fr){ return ap->Z(fr); })
        .def("getXYZ", [](Atom_proxy* ap, int fr){ return ap->XYZ(fr); })
        .def("setX", [](Atom_proxy* ap, int fr, float data){ ap->X(fr) = data; })
        .def("setY", [](Atom_proxy* ap, int fr, float data){ ap->Y(fr) = data; })
        .def("setZ", [](Atom_proxy* ap, int fr, float data){ ap->Z(fr) = data; })
        .def("setXYZ", [](Atom_proxy* ap, int fr, Vector3f_const_ref data){ ap->XYZ(fr) = data; })
    ;

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

        // Get and set
        .def("get_frame",&Selection::get_frame)
        .def("get_system",&Selection::get_system, py::return_value_policy::reference_internal)
        .def("get_text",&Selection::get_text)
        .def("get_index",&Selection::get_index)

        .def("get_chain",&Selection::get_chain)
        .def("get_unique_chain",&Selection::get_unique_chain)
        .def("set_chain",py::overload_cast<char>(&Selection::set_chain))
        .def("set_chain",py::overload_cast<const std::vector<char>&>(&Selection::set_chain))

        .def("get_resid",&Selection::get_resid)
        .def("get_unique_resid",&Selection::get_unique_resid)
        .def("set_resid",py::overload_cast<int>(&Selection::set_resid))
        .def("set_resid",py::overload_cast<const std::vector<int>&>(&Selection::set_resid))

        .def("get_resindex",&Selection::get_resindex)
        .def("get_unique_resindex",&Selection::get_unique_resindex)

        .def("get_name",&Selection::get_name)
        .def("set_name",py::overload_cast<string&>(&Selection::set_name))
        .def("set_name",py::overload_cast<const std::vector<string>&>(&Selection::set_name))

        .def("get_resname",&Selection::get_resname)
        .def("get_unique_resname",&Selection::get_unique_resname)
        .def("set_resname",py::overload_cast<string&>(&Selection::set_resname))
        .def("set_resname",py::overload_cast<const std::vector<string>&>(&Selection::set_resname))

        .def("get_xyz", [](Selection* sel){ return sel->get_xyz(true); }) // pass true for row-major matrix
        .def("set_xyz", &Selection::set_xyz) // detects raw-major matrix internally

        .def("get_mass",&Selection::get_mass)
        .def("set_mass",py::overload_cast<float>(&Selection::set_mass))
        .def("set_mass",py::overload_cast<const std::vector<float>>(&Selection::set_mass))

        .def("get_beta",&Selection::get_beta)
        .def("set_beta",py::overload_cast<float>(&Selection::set_beta))
        .def("set_beta",py::overload_cast<const std::vector<float>&>(&Selection::set_beta))

        .def("get_occupancy",&Selection::get_occupancy)
        .def("set_occupancy",py::overload_cast<float>(&Selection::set_occupancy))
        .def("set_occupancy",py::overload_cast<const std::vector<float>&>(&Selection::set_occupancy))

        .def("get_tag",&Selection::get_tag)
        .def("set_tag",py::overload_cast<string&>(&Selection::set_tag))
        .def("set_tag",py::overload_cast<const std::vector<string>&>(&Selection::set_tag))

        // Properties
        .def("center",&Selection::center,"mass_weighted"_a=false,"periodic"_a=false)
        .def("minmax",[](Selection* sel){Vector3f min,max; sel->minmax(min,max); return py::make_tuple(min,max);})

        .def("sasa", [](Selection* sel, float probe_r, bool do_total_volume, bool do_area_per_atom, bool do_vol_per_atom){
            float vol;
            std::vector<float> area_per_atom;
            std::vector<float> volume_per_atom;
            float* vol_ptr;
            std::vector<float> *area_per_atom_ptr, *volume_per_atom_ptr;
            vol_ptr = do_total_volume ? &vol : nullptr;
            area_per_atom_ptr = do_area_per_atom ? &area_per_atom : nullptr;
            volume_per_atom_ptr = do_vol_per_atom ? &volume_per_atom : nullptr;
            float a = sel->sasa(probe_r,vol_ptr,area_per_atom_ptr,volume_per_atom_ptr);
            py::list ret;
            ret.append(a);
            if(do_total_volume) ret.append(vol);
            if(do_area_per_atom) ret.append(area_per_atom);
            if(do_vol_per_atom) ret.append(volume_per_atom);
            return ret;
         }, "probe_r"_a=0.14, "do_total_volume"_a=false, "do_area_per_atom"_a=false, "do_vol_per_atom"_a=false)

        .def("average_structure", [](Selection* sel, int b, int e){
                return sel->average_structure(b,e,true); // pass true for row-major matrix
            }, "b"_a=0, "e"_a=-1)

        .def("atom_traj", [](Selection* sel, int i, int b, int e){
                return sel->atom_traj(i,b,e,true); // pass true for row-major matrix
            }, "i"_a, "b"_a=0, "e"_a=-1)

        .def("inertia",[](Selection* sel, bool is_periodic){
                Vector3f m;
                Matrix3f ax;
                sel->inertia(m,ax,is_periodic);
                return py::make_tuple(m,ax.transpose());
            },"is_periodic"_a=false)

        .def("gyration",&Selection::gyration, "periodic"_a=false)

        .def("distance", &Selection::distance, "i"_a, "j"_a, "periodic"_a=true, "dims"_a=Eigen::Vector3i::Ones())
        .def("angle", &Selection::angle, "i"_a, "j"_a, "k"_a, "periodic"_a=true, "dims"_a=Eigen::Vector3i::Ones())
        .def("dihedral", &Selection::dihedral, "i"_a, "j"_a, "k"_a, "l"_a, "periodic"_a=true, "dims"_a=Eigen::Vector3i::Ones())

        // Geometry transforms
        .def("translate", &Selection::translate)
        .def("rotate",py::overload_cast<int,float>(&Selection::rotate))
        .def("rotate",py::overload_cast<int,float,Vector3f_const_ref>(&Selection::rotate))
        .def("rotate",py::overload_cast<Vector3f_const_ref,float,Vector3f_const_ref>(&Selection::rotate))
        .def("rotate",[](Selection* sel, Matrix3f_const_ref m){ sel->rotate(m.transpose()); })
        .def("rotate",py::overload_cast<Vector3f_const_ref,Vector3f_const_ref>(&Selection::rotate))
        .def("wrap", &Selection::wrap, "dims"_a=Eigen::Vector3i::Ones())
        .def("unwrap", &Selection::unwrap, "dims"_a=Eigen::Vector3i::Ones())
        .def("unwrap_bonds", &Selection::unwrap_bonds, "d"_a=0.2, "lead_ind"_a=0, "dims"_a=Eigen::Vector3i::Ones())
        .def("principal_transform", [](Selection* sel, bool periodic){
                Matrix4f m = sel->principal_transform(periodic).matrix().transpose();
                return m;
            }, "periodic"_a=false)

        .def("principal_orient",&Selection::principal_orient, "periodic"_a=false)

        // Fitting and rmsd
        .def("rmsd",py::overload_cast<int>(&Selection::rmsd,py::const_))
        .def("rmsd",py::overload_cast<int,int>(&Selection::rmsd,py::const_))
        .def("fit_trajectory",&Selection::fit_trajectory, "ref_frame"_a=0, "b"_a=0, "e"_a=-1)

        .def("fir_transform", [](Selection* sel, int fr1, int fr2){
                Matrix4f m = sel->fit_transform(fr1,fr2).matrix().transpose();
                return m;
            })

        .def("apply_transform", [](Selection* sel, const Eigen::Ref<const Eigen::Matrix4f>& m){
                Affine3f t(m.transpose());
                sel->apply_transform(t);
            })

        // Energy
        .def("non_bond_energy", &Selection::non_bond_energy, "cutoff"_a=0.25, "periodic"_a=true)

        // IO
        .def("write", py::overload_cast<string,int,int>(&Selection::write), "fname"_a, "b"_a=0, "e"_a=-1)

        // Util
        .def("size",&Selection::size)
        .def("__len__", &Selection::size)
        .def("text_based",&Selection::text_based)
        .def("coord_dependent",&Selection::coord_dependent)
        .def("flatten",&Selection::flatten)
        .def("gromacs_ndx",&Selection::gromacs_ndx)

        // Indexing and iterating
        .def("__iter__", [](Selection* s) {
            return py::make_iterator(s->begin(), s->end());
        }, py::keep_alive<0,1>() /* Essential: keep object alive while iterator exists */)

        .def("__getitem__", [](const Selection &s, size_t i) {
                if (i >= s.size()) throw py::index_error();
                return s[i]; // Returns atom proxy object
            })

        // Splitting
        .def("split_by_connectivity", [](Selection* sel,float d,bool periodic){
                std::vector<Selection> res;
                sel->split_by_connectivity(d,res,periodic);
                return res;
            })

        .def("split_by_residue", [](Selection* sel){
                std::vector<Selection> res;
                sel->split_by_residue(res);
                return res;
            })

        .def("split_by_chain", [](Selection* sel){
                std::vector<Selection> res;
                sel->split_by_chain(res);
                return res;
            })

        .def("split_by_contiguous_index", [](Selection* sel){
                std::vector<Selection> res;
                sel->split_by_contiguous_index(res);
                return res;
            })

        .def("split_by_contiguous_residue", [](Selection* sel){
                std::vector<Selection> res;
                sel->split_by_contiguous_residue(res);
                return res;
            })

        .def("each_residue", [](Selection* sel){
                std::vector<Selection> res;
                sel->each_residue(res);
                return res;
            })

        // split based on callback have to be implemented in python side
        // since no means to bind templated return value!

        // Accessors
        .def("VDW",&Selection::VDW)
        .def_property("box", [](Selection* obj){return obj->Box();}, [](Selection* obj,const Periodic_box& val){obj->Box()=val;})
        .def_property("time", [](Selection* obj){return obj->Time();}, [](Selection* obj, float val){obj->Time()=val;})

        // No other accessors are exposed in favor to [] operator
    ;

    // Free functions
    m.def("rmsd",[](const Selection& sel1, const Selection& sel2){ return rmsd(sel1,sel2); });
    m.def("rmsd",[](const Selection& sel1, int fr1, const Selection& sel2, int fr2){ return rmsd(sel1,fr1,sel2,fr2); });
    m.def("fit",[](Selection& sel1, const Selection& sel2){ fit(sel1,sel2); });
    m.def("non_bond_energy", [](const Selection& sel1, const Selection& sel2,float cutoff,int fr,bool periodic){
        return non_bond_energy(sel1,sel2,cutoff,fr,periodic);
    },"sel1"_a, "sel2"_a, "cutoff"_a=0.25, "fr"_a=-1, "periodic"_a=true);
    m.def("copy_coord",[](const Selection& sel1, int fr1, Selection& sel2, int fr2){ return copy_coord(sel1,fr1,sel2,fr2); });

}
