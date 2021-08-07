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
#include "pteros/core/pteros_error.h"
#include "bindings_util.h"

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace Eigen;
using namespace pybind11::literals;

#define DEF_PROPERTY(_name,_dtype) \
    .def_property(#_name, [](AtomHandler* obj){return obj->_name();}, [](AtomHandler* obj,const _dtype& val){obj->_name()=val;})

void make_bindings_Selection(py::module& m){

    using RowMatrixXf = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    py::class_<AtomHandler>(m, "AtomHandler")
        DEF_PROPERTY(resid,int)
        DEF_PROPERTY(resindex,int)
        DEF_PROPERTY(resname,string)
        DEF_PROPERTY(name,string)
        DEF_PROPERTY(chain,char)
        DEF_PROPERTY(tag,string)
        DEF_PROPERTY(occupancy,float)
        DEF_PROPERTY(beta,float)
        DEF_PROPERTY(mass,float)
        DEF_PROPERTY(charge,float)
        DEF_PROPERTY(type,int)
        DEF_PROPERTY(atomic_number,int)
        DEF_PROPERTY(type_name,string)
        DEF_PROPERTY(atom,Atom)
        DEF_PROPERTY(x,float)
        DEF_PROPERTY(y,float)
        DEF_PROPERTY(z,float)
        .def_property("xyz", [](AtomHandler* obj){return obj->xyz();}, [](AtomHandler* obj,Vector3f_const_ref val){obj->xyz()=val;})
        //.def_property("vel", [](AtomHandler* obj){return obj->vel();}, [](AtomHandler* obj,Vector3f_const_ref val){obj->vel()=val;})
        //.def_property("force", [](AtomHandler* obj){return obj->force();}, [](AtomHandler* obj,Vector3f_const_ref val){obj->force()=val;})
        .def_property_readonly("element_name", [](AtomHandler* obj){return obj->element_name();})
        .def_property_readonly("vdw", [](AtomHandler* obj){return obj->vdw();})        
        .def_property_readonly("index", [](AtomHandler* obj){return obj->index();})
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
        .def("__call__", py::overload_cast<string>(&Selection::operator()), py::keep_alive<0,1>())
        .def("__call__", py::overload_cast<int,int>(&Selection::operator()), py::keep_alive<0,1>())
        .def("__call__", py::overload_cast<const std::vector<int>&>(&Selection::operator()), py::keep_alive<0,1>())
        .def("select", py::overload_cast<string>(&Selection::operator()), py::keep_alive<0,1>())
        .def("select", py::overload_cast<int,int>(&Selection::operator()), py::keep_alive<0,1>())
        .def("select", py::overload_cast<const std::vector<int>&>(&Selection::operator()), py::keep_alive<0,1>())

        // Get and set
        .def("get_frame",&Selection::get_frame)
        .def("set_frame",&Selection::set_frame)
        .def("get_system",&Selection::get_system, py::return_value_policy::reference_internal)
        .def("get_text",&Selection::get_text)
        .def("get_index",&Selection::get_index)        

        .def("get_chain",&Selection::get_chain,"unique"_a=false)
        .def("set_chain",py::overload_cast<char>(&Selection::set_chain))
        .def("set_chain",py::overload_cast<const std::vector<char>&>(&Selection::set_chain))

        .def("get_resid",&Selection::get_resid,"unique"_a=false)
        .def("set_resid",py::overload_cast<int>(&Selection::set_resid))
        .def("set_resid",py::overload_cast<const std::vector<int>&>(&Selection::set_resid))

        .def("get_resindex",&Selection::get_resindex,"unique"_a=false)

        .def("get_name",&Selection::get_name,"unique"_a=false)
        .def("set_name",py::overload_cast<string>(&Selection::set_name))
        .def("set_name",py::overload_cast<const std::vector<string>&>(&Selection::set_name))

        .def("get_resname",&Selection::get_resname,"unique"_a=false)
        .def("set_resname",py::overload_cast<string>(&Selection::set_resname))
        .def("set_resname",py::overload_cast<const std::vector<string>&>(&Selection::set_resname))

        .def("get_xyz", [](Selection* sel){ return sel->get_xyz(true); }) // pass true for row-major matrix
        .def("set_xyz", &Selection::set_xyz) // detects raw-major matrix internally

        .def("get_vel", [](Selection* sel){ return sel->get_vel(true); }) // pass true for row-major matrix
        .def("set_vel", &Selection::set_vel) // detects raw-major matrix internally

        .def("get_force", [](Selection* sel){ return sel->get_force(true); }) // pass true for row-major matrix
        .def("set_force", &Selection::set_force) // detects raw-major matrix internally

        .def("get_mass",&Selection::get_mass)
        .def("set_mass",py::overload_cast<float>(&Selection::set_mass))
        .def("set_mass",py::overload_cast<const std::vector<float>&>(&Selection::set_mass))

        .def("get_beta",&Selection::get_beta)
        .def("set_beta",py::overload_cast<float>(&Selection::set_beta))
        .def("set_beta",py::overload_cast<const std::vector<float>&>(&Selection::set_beta))

        .def("get_occupancy",&Selection::get_occupancy)
        .def("set_occupancy",py::overload_cast<float>(&Selection::set_occupancy))
        .def("set_occupancy",py::overload_cast<const std::vector<float>&>(&Selection::set_occupancy))

        .def("get_charge",&Selection::get_charge)
        .def("get_total_charge",&Selection::get_total_charge)
        .def("set_charge",py::overload_cast<float>(&Selection::set_charge))
        .def("set_charge",py::overload_cast<const std::vector<float>&>(&Selection::set_charge))

        .def("get_tag",&Selection::get_tag,"unique"_a=false)
        .def("set_tag",py::overload_cast<string>(&Selection::set_tag))
        .def("set_tag",py::overload_cast<const std::vector<string>&>(&Selection::set_tag))

        // Properties
        .def("center",py::overload_cast<bool,Array3i_const_ref,int>(&Selection::center,py::const_),"mass_weighted"_a=false,"pbc"_a=noPBC,"pbc_atom"_a=-1)
        .def("center",py::overload_cast<const std::vector<float>&,Array3i_const_ref,int>(&Selection::center,py::const_),"weight"_a,"pbc"_a=noPBC,"pbc_atom"_a=-1)

        .def("minmax",[](Selection* sel){Vector3f min,max; sel->minmax(min,max); return py::make_tuple(min,max);})
        .def("is_large",&Selection::is_large)

        .def("powersasa", [](Selection* sel, float probe_r, bool do_area_per_atom, bool do_total_volume, bool do_vol_per_atom){
            float vol;
            std::vector<float> area_per_atom;
            std::vector<float> volume_per_atom;
            float* vol_ptr;
            std::vector<float> *area_per_atom_ptr, *volume_per_atom_ptr;
            vol_ptr = do_total_volume ? &vol : nullptr;
            area_per_atom_ptr = do_area_per_atom ? &area_per_atom : nullptr;
            volume_per_atom_ptr = do_vol_per_atom ? &volume_per_atom : nullptr;
            float a = sel->powersasa(probe_r,area_per_atom_ptr,vol_ptr,volume_per_atom_ptr);
            py::list ret;
            ret.append(a);            
            if(do_area_per_atom) ret.append(area_per_atom);
            if(do_total_volume) ret.append(vol);
            if(do_vol_per_atom) ret.append(volume_per_atom);
            return ret;
         }, "probe_r"_a=0.14, "do_area_per_atom"_a=false, "do_total_volume"_a=false, "do_vol_per_atom"_a=false)

        .def("sasa", [](Selection* sel, float probe_r, bool do_area_per_atom, int n_sphere_points){
            std::vector<float> area_per_atom;
            std::vector<float> *area_per_atom_ptr;
            area_per_atom_ptr = do_area_per_atom ? &area_per_atom : nullptr;
            float a = sel->sasa(probe_r,area_per_atom_ptr,n_sphere_points);
            py::list ret;
            ret.append(a);
            if(do_area_per_atom) ret.append(area_per_atom);
            return ret;
        }, "probe_r"_a=0.14, "do_area_per_atom"_a=false, "n_sphere_points"_a=960)

        .def("average_structure", [](Selection* sel, int b, int e){
                return sel->average_structure(b,e,true); // pass true for row-major matrix
            }, "b"_a=0, "e"_a=-1)

        .def("atom_traj", [](Selection* sel, int i, int b, int e){
                return sel->atom_traj(i,b,e,true); // pass true for row-major matrix
            }, "i"_a, "b"_a=0, "e"_a=-1)

        .def("inertia",[](Selection* sel, Array3i_const_ref pbc, bool pbc_atom){
                Vector3f m;
                Matrix3f ax;
                sel->inertia(m,ax,pbc,pbc_atom);
                return py::make_tuple(m,ax.transpose());
            },"pbc"_a=noPBC,"pbc_atom"_a=-1)

        .def("gyration",&Selection::gyration, "pbc"_a=noPBC,"pbc_atom"_a=-1)
        .def("dipole",&Selection::dipole, "is_charged"_a=false,"pbc"_a=noPBC,"pbc_atom"_a=-1)

        .def("distance", &Selection::distance, "i"_a, "j"_a, "pbc"_a=fullPBC)
        .def("angle", &Selection::angle, "i"_a, "j"_a, "k"_a, "pbc"_a=fullPBC)
        .def("dihedral", &Selection::dihedral, "i"_a, "j"_a, "k"_a, "l"_a, "pbc"_a=fullPBC)
        .def("num_residues",&Selection::num_residues)

        // Geometry transforms
        .def("translate", &Selection::translate)
        .def("translate_to", &Selection::translate_to, "vec"_a, "mass_weighted"_a=false, "pbc"_a=noPBC, "pbc_atom"_a=-1)
        .def("rotate",&Selection::rotate, "pivot"_a, "axis"_a, "angle"_a)

        .def("wrap", &Selection::wrap, "pbc"_a=fullPBC)
        .def("unwrap", &Selection::unwrap, "pbc"_a=fullPBC, "pbc_atom"_a=-1)
        .def("unwrap_bonds", &Selection::unwrap_bonds, "d"_a, "pbc"_a=fullPBC, "pbc_atom"_a=-1)
        .def("principal_transform", [](Selection* sel, Array3i_const_ref pbc, bool pbc_atom){
                Matrix4f m = sel->principal_transform(pbc,pbc_atom).matrix().transpose();
                return m;
            }, "pbc"_a=noPBC,"pbc_atom"_a=-1)

        .def("principal_orient",&Selection::principal_orient,"pbc"_a=noPBC,"pbc_atom"_a=-1)

        // Fitting and rmsd
        .def("rmsd",py::overload_cast<int>(&Selection::rmsd,py::const_))
        .def("rmsd",py::overload_cast<int,int>(&Selection::rmsd,py::const_))
        .def("fit_trajectory",&Selection::fit_trajectory, "ref_frame"_a=0, "b"_a=0, "e"_a=-1)
        .def("fit",&Selection::fit)

        .def("fit_transform", [](Selection* sel, int fr1, int fr2){
                Matrix4f m = sel->fit_transform(fr1,fr2).matrix().transpose();
                return m;
            })

        .def("apply_transform", [](Selection* sel, const Eigen::Ref<const Eigen::Matrix4f>& m){
                Affine3f t(m.transpose());
                sel->apply_transform(t);
            })

        // Energy
        .def("non_bond_energy", &Selection::non_bond_energy, "cutoff"_a=0.0, "pbc"_a=true)

        // IO
        .def("write", py::overload_cast<string,int,int>(&Selection::write,py::const_), "fname"_a, "b"_a=-1, "e"_a=-1)

        // Util
        .def("is_large",&Selection::is_large)
        .def("size",&Selection::size)
        .def("__len__", &Selection::size)
        .def("text_based",&Selection::text_based)
        .def("coord_dependent",&Selection::coord_dependent)
        .def("flatten",&Selection::flatten)
        .def("to_gromacs_ndx",&Selection::to_gromacs_ndx)
        .def("find_index",&Selection::find_index)

        // Indexing and iterating
        .def("__iter__", [](Selection* s) {
            return py::make_iterator(s->begin(), s->end());
        }, py::keep_alive<0,1>() /* Essential: keep object alive while iterator exists */)

        .def("__getitem__", [](Selection &s, size_t i) {
                if(i >= s.size()) throw py::index_error();                
                return s[i]; // Returns atom proxy object
            }, py::keep_alive<0,1>())

        .def("__getitem__", [](Selection &s, py::tuple ind_fr) {
                int i = ind_fr[0].cast<int>();
                int fr = ind_fr[1].cast<int>();
                if(i >= s.size() || fr<0 || fr>=s.get_system()->num_frames()) throw py::index_error();
                return s[{i,fr}]; // Returns atom proxy object
            }, py::keep_alive<0,1>())

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

        .def("split_by_molecule", [](Selection* sel){
                std::vector<Selection> res;
                sel->split_by_molecule(res);
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

        // split based on callback have to be implemented on python side
        // since no means to bind templated return value!

        // dssp
        .def("dssp", py::overload_cast<string>(&Selection::dssp, py::const_))
        .def("dssp", py::overload_cast<>(&Selection::dssp, py::const_))

        // Accessors
        .def_property("box", [](Selection* obj){return obj->box();}, [](Selection* obj,const PeriodicBox& val){obj->box()=val;})
        .def_property("time", [](Selection* obj){return obj->time();}, [](Selection* obj, float val){obj->time()=val;})

        // No other accessors are exposed in favor to [] operator
    ;

    // Free functions
    m.def("rmsd",[](const Selection& sel1, const Selection& sel2){ return rmsd(sel1,sel2); });
    m.def("rmsd",[](const Selection& sel1, int fr1, const Selection& sel2, int fr2){ return rmsd(sel1,fr1,sel2,fr2); });

    m.def("fit",[](Selection& sel1, const Selection& sel2){ fit(sel1,sel2); });
    m.def("fit_transform",[](Selection& sel1, const Selection& sel2){
        Matrix4f m = fit_transform(sel1,sel2).matrix().transpose();
        return m;
    });

    m.def("non_bond_energy", [](const Selection& sel1, const Selection& sel2,float cutoff,int fr,bool pbc){
        return non_bond_energy(sel1,sel2,cutoff,fr,pbc);
    },"sel1"_a, "sel2"_a, "cutoff"_a=0.0, "fr"_a=-1, "pbc"_a=true);

    m.def("copy_coord",[](const Selection& sel1, int fr1, Selection& sel2, int fr2){ return copy_coord(sel1,fr1,sel2,fr2); });
    m.def("copy_coord",[](const Selection& sel1, Selection& sel2){ return copy_coord(sel1,sel2); });
}




