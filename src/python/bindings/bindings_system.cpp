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

#include "bindings_system.h"
#include "pteros/core/system.h"
#include "pteros/python/bindings_util.h"

using namespace pteros;
using namespace Eigen;
namespace bp = boost::python;
using namespace boost::python;

/**********************
  Wrappers for System
***********************/

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(load_overloads, load, 1, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(frame_delete_overloads, frame_delete, 0, 2)

const Periodic_box& System_getBox(System* s, int fr){
    return s->Box(fr);
}

void System_setBox(System* s, int fr, Periodic_box& b){
     s->Box(fr) = b;
}

float System_getTime(System* s, int fr){
    return s->Time(fr);
}

void System_setTime(System* s, int fr, float t){
    s->Time(fr) = t;
}

PyObject* System_getXYZ(System* s, int ind, int fr){
    CREATE_PYARRAY_1D_AND_MAP(p,Vector3f,v,3)
    v = s->XYZ(ind,fr);
    return boost::python::incref(p);
}

void System_setXYZ(System* s, PyObject* arr, int ind, int fr){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,v,arr)
     s->XYZ(ind,fr) = v;
}

Frame System_getFrame_data(System* s, int fr){
    return s->Frame_data(fr);
}

void System_setFrame_data(System* s, Frame& data, int fr){
    s->Frame_data(fr) = data;
}


Selection System_atoms_dup(System* s, boost::python::list& data){
    vector<int> r;
    r.resize(len(data));
    for(int i=0;i<r.size();++i) r[i] = extract<int>(data[i]);
    return s->atoms_dup(r);
}

Selection System_atoms_add(System* s, boost::python::list& atm,
                       boost::python::list& crd){
    vector<Atom> a;
    vector<Vector3f> c;
    a.resize(len(atm));
    c.resize(len(crd));
    for(int i=0;i<a.size();++i) a[i] = extract<Atom>(atm[i]);
    for(int i=0;i<c.size();++i){
        boost::python::object o = crd[i];
        MAP_EIGEN_TO_PYTHON_F(Vector3f,v,o.ptr())
        c[i] = v;
    }
    return s->atoms_add(a,c);
}

Selection System_append(System* s, const Atom& atm, PyObject* crd){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,coord,crd)
    return s->append(atm,coord);
}


void System_load_callback(System* sys, string fname, int b, int e, int skip, boost::python::object obj){
    // Create a callback from obj
    auto callback = [&obj](System* s, int fr)->bool { return extract<bool>(obj(ptr(s),fr)); };
    sys->load(fname,b,e,skip,callback);
}

void System_wrap_all1(System* sys, int fr, bp::list& dims){
    Vector3i d;
    for(int i=0;i<3;++i) d(i) = extract<int>(dims[i]);
    sys->wrap_all(fr,d);
}

void System_wrap_all2(System* sys, int fr){
    sys->wrap_all(fr);
}


Selection System_select_vec(System* sys, bp::list l){
    vector<int> v(len(l));
    for(int i=0;i<len(l);++i) v[i] = extract<int>(l[i]);
    return sys->select(v);
}

void make_bindings_System(){
    import_array();

    class_<System, boost::noncopyable>("System", init<>())
        .def(init<std::string>() )
        .def(init<System>() )            
        .def("num_atoms", &System::num_atoms)
        .def("num_frames", &System::num_frames)

        .def("select", static_cast<Selection(System::*) (std::string)>         (&System::select))
        .def("select", static_cast<Selection(System::*) (int ind1, int ind2)>  (&System::select))
        .def("select", &System_select_vec)
        .def("select_all", &System::select_all)

        .def("load", &System::load, load_overloads())
        .def("load", &System_load_callback,(bp::arg("fname"),bp::arg("b")=0,bp::arg("e")=-1,bp::arg("skip")=0,bp::arg("on_frame")))

        .def("frame_append", &System::frame_append)
        .def("frame_dup", &System::frame_dup)
        .def("frame_copy", &System::frame_copy)
        .def("frame_delete", &System::frame_delete, frame_delete_overloads())

        .def("getFrame_data", &System_getFrame_data)
        .def("setFrame_data", &System_setFrame_data)        

        .def("getBox", &System_getBox,return_value_policy<reference_existing_object>())
        .def("setBox", &System_setBox)        

        .def("getTime", &System_getTime)
        .def("setTime", &System_setTime)               

        .def("getXYZ", &System_getXYZ)
        .def("setXYZ", &System_setXYZ)        

        .def("assign_resindex", &System::assign_resindex)

        .def("atoms_dup", &System_atoms_dup)
        .def("atoms_add", &System_atoms_add)

        .def("wrap_all", &System_wrap_all1)
        .def("wrap_all", &System_wrap_all2)

        .def("append", static_cast<Selection(System::*)(const Selection&)>(&System::append))
        .def("append", static_cast<Selection(System::*)(const System&)>(&System::append))
        .def("append", &System_append)

        .def("dssp", static_cast<void(System::*)(std::string)const>(&System::dssp))
        .def("dssp", static_cast<std::string(System::*)()const>(&System::dssp))

        .def("sort_by_resindex",&System::sort_by_resindex)
    ;
}
