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
using namespace std;
namespace bp = boost::python;
using namespace boost::python;

/**********************
  Wrappers for System
***********************/

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
    CREATE_PYARRAY_1D_AND_MAP_F(p,Vector3f,v,3)
    v = s->XYZ(ind,fr);
    return p;
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

Selection System_atoms_delete(System* s, boost::python::list& atm){
    vector<int> a;
    a.resize(len(atm));
    for(int i=0;i<a.size();++i) a[i] = extract<int>(atm[i]);
    s->atoms_delete(a);
}


void System_rearrange(System* s, boost::python::list& data){
    // Try to extract string from the first element
    if(len(data)==0) throw Pteros_error("Need a list of selection or the list of selection strings!");
    extract<string> get_str(data[0]);
    if (get_str.check()){
        // Work with strings
        vector<string> sels(len(data));
        for(int i=0;i<sels.size();++i) sels[i] = extract<string>(data[i]);
        s->rearrange(sels);
    } else {
        // Work with selections
        vector<Selection> sels(len(data));
        for(int i=0;i<sels.size();++i) sels[i] = extract<Selection>(data[i]);
        s->rearrange(sels);
    }
}

void System_load_normal(System* sys, string fname, int b=0, int e=-1,int skip=0){
    sys->load(fname,b,e,skip);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(system_load_overloads, System_load_normal, 2, 5)

void System_load_callback(System* sys, string fname, int b, int e, int skip, boost::python::object obj){
    // Create a callback from obj
    auto callback = [&obj](System* s, int fr)->bool { return extract<bool>(obj(ptr(s),fr)); };
    sys->load(fname,b,e,skip,callback);
}

void System_wrap_all1(System* sys, int fr, bp::list& dims){
    Vector3i d;
    for(int i=0;i<3;++i) d(i) = extract<int>(dims[i]);
    sys->wrap(fr,d);
}

void System_wrap_all2(System* sys, int fr){
    sys->wrap(fr);
}

float System_distance1(System* sys, int i, int j, int fr, bool pbc, PyObject* dims){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims)
    return sys->distance(i,j,fr,pbc,dim);
}

float System_distance2(System* sys, int i, int j,  int fr, bool pbc){
    return sys->distance(i,j,fr,pbc);
}

float System_distance3(System* sys, int i, int j, int fr){
    return sys->distance(i,j,fr);
}

float System_angle1(System* sys, int i, int j, int k,  int fr, bool pbc, PyObject* dims){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims)
    return sys->angle(i,j,k,fr,pbc,dim);
}

float System_angle2(System* sys, int i, int j, int k,  int fr, bool pbc){
    return sys->angle(i,j,k,fr,pbc);
}

float System_angle3(System* sys, int i, int j, int k, int fr){
    return sys->angle(i,j,k,fr);
}

float System_dihedral1(System* sys, int i, int j, int k, int l, int fr, bool pbc, PyObject* dims){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims)
    return sys->dihedral(i,j,k,l,fr,pbc,dim);
}

float System_dihedral2(System* sys, int i, int j, int k, int l, int fr, bool pbc){
    return sys->dihedral(i,j,k,l,fr,pbc);
}

float System_dihedral3(System* sys, int i, int j, int k, int l, int fr){
    return sys->dihedral(i,j,k,l,fr);
}

void System_distribute(System* s, const Selection& sel, boost::python::list& ncopy, PyObject* shift){
    MAP_EIGEN_TO_PYTHON_F(Matrix3f,sh,shift)
    Vector3i nc;
    for(int i=0;i<3;++i) nc(i) = extract<int>(ncopy[i]);
    s->distribute(sel,nc,sh);
}

void set_filter_from_py_obj(System* s, PyObject* obj){
    if( PyString_Check(obj) ) {
        string str = extract<string>(object( handle<>(borrowed(obj)) ));
        s->set_filter(str);
    } else if( PySequence_Check(obj) ) {
        bp::object l( handle<>(borrowed(obj)) );
        vector<int> v(len(l));
        for(int i=0;i<len(l);++i) v[i] = extract<int>(l[i]);
        s->set_filter(v);
    } else {
        throw Pteros_error("Invalid arguments for set_filter!");
    }
}


Selection system_select_2_args(System* sys, PyObject* obj, int fr){
    if(PyCallable_Check(obj)){
        // Create a callback from obj
        auto callback = [&obj](const System& s, int fr, vector<int>& ind)->void {
            auto l = bp::call<bp::list>(obj,boost::ref(s),fr);
            ind.resize(len(l));
            for(int i=0;i<len(l);++i) ind[i] = extract<int>(l[i]);
        };
        return sys->select(callback,fr);
    } else if(PyString_Check(obj)) {
        string s = bp::extract<string>( object(handle<>(borrowed(obj))) );
        return sys->select(s,fr);
    } else if(PyInt_Check(obj)) {
        int i = bp::extract<int>( object(handle<>(borrowed(obj))) );
        return sys->select(i,fr);
    }
}

Selection system_select_1_arg(System* sys, PyObject* obj){
    if(PyCallable_Check(obj)){
        // Create a callback from obj
        auto callback = [&obj](const System& s, int fr, vector<int>& ind)->void {
            auto l = bp::call<bp::list>(obj,boost::ref(s),fr);
            ind.resize(len(l));
            for(int i=0;i<len(l);++i) ind[i] = extract<int>(l[i]);
        };
        return sys->select(callback);
    } else if(PyString_Check(obj)) {
        string s = bp::extract<string>( object(handle<>(borrowed(obj))) );
        return sys->select(s);
    } else if(PySequence_Check(obj)) {
        bp::object l( handle<>(borrowed(obj)) );
        vector<int> ind(len(l));
        for(int i=0;i<len(l);++i) ind[i] = extract<int>(l[i]);
        return sys->select(ind);
    }
}


//==================================================================

void make_bindings_System(){
    import_array();

    class_<System, boost::noncopyable>("System", init<>())
        .def(init<std::string>() )
        .def(init<System>() )            
        .def("num_atoms", &System::num_atoms)
        .def("num_frames", &System::num_frames)

        .def("select",&system_select_1_arg,with_custodian_and_ward_postcall<0,1>())
        .def("__call__",&system_select_1_arg,with_custodian_and_ward_postcall<0,1>())

        .def("select",&system_select_2_args,with_custodian_and_ward_postcall<0,1>())
        .def("__call__",&system_select_2_args,with_custodian_and_ward_postcall<0,1>())

        .def("select_all",&System::select_all,with_custodian_and_ward_postcall<0,1>())
        .def("__call__",&System::select_all,with_custodian_and_ward_postcall<0,1>())

        .def("load", &System_load_normal, system_load_overloads())
        .def("load", &System_load_callback,(bp::arg("fname"),bp::arg("b")=0,bp::arg("e")=-1,bp::arg("skip")=0,bp::arg("on_frame")))

        .def("frame_append", &System::frame_append)
        .def("frame_dup", &System::frame_dup)
        .def("frame_copy", &System::frame_copy)
        .def("frame_swap", &System::frame_swap)
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
        .def("atoms_delete", &System_atoms_delete)

        .def("wrap", &System_wrap_all1)
        .def("wrap", &System_wrap_all2)

        .def("append", static_cast<Selection(System::*)(const Selection&)>(&System::append))
        .def("append", static_cast<Selection(System::*)(const System&)>(&System::append))
        .def("append", &System_append)
        .def("rearrange", &System_rearrange)

        .def("keep", static_cast<void(System::*)(const Selection&)>(&System::keep))
        .def("keep", static_cast<void(System::*)(const string&)>(&System::keep))
        .def("remove", static_cast<void(System::*)(Selection&)>(&System::remove))
        .def("remove", static_cast<void(System::*)(const string&)>(&System::remove))
        .def("distribute",&System_distribute)

        .def("dssp", static_cast<void(System::*)(std::string,int)const>(&System::dssp))
        .def("dssp", static_cast<std::string(System::*)(int)const>(&System::dssp))

        .def("sort_by_resindex",&System::sort_by_resindex)

        .def("distance",&System_distance1)
        .def("distance",&System_distance2)
        .def("distance",&System_distance3)

        .def("angle",&System_angle1)
        .def("angle",&System_angle2)
        .def("angle",&System_angle3)

        .def("dihedral",&System_dihedral1)
        .def("dihedral",&System_dihedral2)
        .def("dihedral",&System_dihedral3)

        .def("set_filter", static_cast<void(System::*)(int,int)>(&System::set_filter))
        .def("set_filter", &set_filter_from_py_obj)

        .def("clear",&System::clear)
    ;
}
