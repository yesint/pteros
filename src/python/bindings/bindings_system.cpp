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





#include "pteros/core/selection.h"
#include "bindings_util.h"

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace Eigen;
using namespace pybind11::literals;

void make_bindings_System(py::module& m){

    py::class_<System>(m, "System")
        .def(py::init<>())
        .def(py::init<const std::string &>())
        .def(py::init<const Selection &>())

        // Size
        .def("num_atoms", &System::num_atoms)
        .def("num_frames", &System::num_frames)

        // Append
        .def("append", py::overload_cast<const System&>(&System::append))
        .def("append", py::overload_cast<const Selection&,bool>(&System::append),"sel"_a,"current_frame"_a=false)
        .def("append", py::overload_cast<const Atom&, Vector3f_const_ref>(&System::append))
        .def("append", py::overload_cast<const AtomHandler&>(&System::append))

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

        // Writing
        .def("write", py::overload_cast<string,int,int>(&System::write,py::const_), "fname"_a, "b"_a=0, "e"_a=-1)

        // Selecting
        .def("__call__", py::overload_cast<>(&System::operator()), py::keep_alive<0,1>())
        .def("__call__", py::overload_cast<string,int>(&System::operator()),"str"_a,"fr"_a=0, py::keep_alive<0,1>())
        .def("__call__", py::overload_cast<int,int>(&System::operator()), py::keep_alive<0,1>())
        .def("__call__", py::overload_cast<const std::vector<int>&>(&System::operator()), py::keep_alive<0,1>())
        .def("__call__", py::overload_cast<const std::function<void(const System&,int,std::vector<int>&)>&,int>(&System::operator()),"callback"_a,"fr"_a=0, py::keep_alive<0,1>())
        .def("select_all", py::overload_cast<>(&System::operator()), py::keep_alive<0,1>())
        .def("select", py::overload_cast<string,int>(&System::operator()),"str"_a,"fr"_a=0, py::keep_alive<0,1>())
        .def("select", py::overload_cast<int,int>(&System::operator()), py::keep_alive<0,1>())
        .def("select", py::overload_cast<const std::vector<int>&>(&System::operator()), py::keep_alive<0,1>())
        .def("select", py::overload_cast<const std::function<void(const System&,int,std::vector<int>&)>&,int>(&System::operator()),"callback"_a,"fr"_a=0, py::keep_alive<0,1>())

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

        //.def("getBox", py::overload_cast<int>(&System::box, py::const_), "fr"_a=0)
        .def("getBox", py::overload_cast<int>(&System::box), "fr"_a=0, py::return_value_policy::reference_internal)

        .def("setBox", [](System* s,const PeriodicBox& b, int fr){ s->box(fr)=b; }, "box"_a, "fr"_a=0)

        .def("getTime", py::overload_cast<int>(&System::time, py::const_), "fr"_a=0)
        .def("setTime", [](System* s,int t,int fr){ s->time(fr)=t; }, "t"_a, "fr"_a=0)

        .def("getFrame", py::overload_cast<int>(&System::frame, py::const_))
        .def("setFrame", [](System* s, int i, const Frame& fr){ s->frame(i)=fr; })

        .def("getXYZ", py::overload_cast<int,int>(&System::xyz, py::const_), "i"_a, "fr"_a=0)
        .def("setXYZ", [](System* s,Vector3f_const_ref v,int i,int fr){ s->xyz(i,fr)=v; }, "xyz"_a, "i"_a, "fr"_a=0)

        .def("getForce", py::overload_cast<int,int>(&System::force, py::const_), "i"_a, "fr"_a=0)
        .def("setForce", [](System* s,Vector3f_const_ref v,int i,int fr){ s->force(i,fr)=v; }, "force"_a, "i"_a, "fr"_a=0)

        .def("getVel", py::overload_cast<int,int>(&System::vel, py::const_), "i"_a, "fr"_a=0)
        .def("setVel", [](System* s,Vector3f_const_ref v,int i,int fr){ s->vel(i,fr)=v; }, "vel"_a, "i"_a, "fr"_a=0)

        .def("getAtom", py::overload_cast<int>(&System::atom, py::const_))
        .def("setAtom", [](System* s, int i, const Atom& a){ s->atom(i)=a; })

        // operations with atoms
        .def("atoms_dup", &System::atoms_dup)
        .def("atoms_add", &System::atoms_add)
        .def("atoms_delete", &System::atoms_delete)
        .def("atom_move", &System::atom_move)
        .def("atom_clone", &System::atom_clone)
        .def("atom_add_1h", &System::atom_add_1h,"target"_a,"at1"_a,"at2"_a,"at3"_a,"dist"_a=0.109,"pbc"_a=true)
        .def("atom_add_2h", &System::atom_add_2h,"target"_a,"at1"_a,"at2"_a,"dist"_a=0.109,"pbc"_a=true)
        .def("atom_add_3h", &System::atom_add_3h,"target"_a,"at1"_a,"dist"_a=0.109,"pbc"_a=true)

        // wrap
        .def("wrap", &System::wrap, "fr"_a, "pbc"_a=fullPBC)

        // measuring
        .def("distance", &System::distance, "i"_a, "j"_a, "fr"_a, "pbc"_a=fullPBC)
        .def("angle", &System::angle, "i"_a, "j"_a, "k"_a, "fr"_a, "pbc"_a=fullPBC)
        .def("dihedral", &System::dihedral, "i"_a, "j"_a, "k"_a, "l"_a, "fr"_a, "pbc"_a=fullPBC)

         // Util
        .def("clear", &System::clear)
        .def("clear_vel", &System::clear_vel)
        .def("clear_force", &System::clear_force)
        .def("force_field_ready", &System::force_field_ready)
        .def("assign_resindex", &System::assign_resindex, "start"_a=0)
        .def("sort_by_resindex", &System::sort_by_resindex)

        // Indexing
        .def("__getitem__", [](System &s, py::tuple ind_fr) {
                int i = ind_fr[0].cast<int>();
                int fr = ind_fr[1].cast<int>();
                if(i >= s.num_atoms() || fr<0 || fr>=s.num_frames()) throw py::index_error();
                return s[{i,fr}]; // Returns atom proxy object
            }, py::keep_alive<0,1>())

    ;
}




