/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
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

/**********************
  Wrappers for System
***********************/

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(load_overloads, load, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(frame_delete_overloads, frame_delete, 0, 2)

PyObject* System_getBox(System* s, int fr){        
    CREATE_PYARRAY_2D(p,3,3)
    MAP_EIGEN_TO_PYARRAY(m,Matrix3f,p)    
    m = s->Box(fr);
    return boost::python::incref(p);
}
void System_setBox(System* s, PyObject* arr, int fr){
    MAP_EIGEN_TO_PYARRAY(m,Matrix3f,arr)
     s->Box(fr) = m;
}

void System_get_box_vectors_angles(System* s, int fr, PyObject* vectors, PyObject* angles){
    MAP_EIGEN_TO_PYARRAY(v,Vector3f,vectors)
    MAP_EIGEN_TO_PYARRAY(a,Vector3f,angles)
    Vector3f vv,aa;
    s->get_box_vectors_angles(fr,vv,aa);
    v = vv;
    a = aa;
}

float System_getTime(System* s, int fr){
    return s->Time(fr);
}

void System_setTime(System* s, int fr, float t){
    s->Time(fr) = t;
}

PyObject* System_getXYZ(System* s, int ind, int fr){
    CREATE_PYARRAY_1D(p,3)
    MAP_EIGEN_TO_PYARRAY(v,Vector3f,p)
    v = s->XYZ(ind,fr);
    return boost::python::incref(p);
}

void System_setXYZ(System* s, PyObject* arr, int ind, int fr){
    MAP_EIGEN_TO_PYARRAY(v,Vector3f,arr)
     s->XYZ(ind,fr) = v;
}

Frame System_getFrame_data(System* s, int fr){
    return s->Frame_data(fr);
}

void System_setFrame_data(System* s, Frame& data, int fr){
    s->Frame_data(fr) = data;
}


void System_atoms_dup1(System* s, boost::python::list& data, Selection* sel){
    vector<int> r;
    r.resize(len(data));
    for(int i=0;i<r.size();++i) r[i] = extract<int>(data[i]);
    s->atoms_dup(r,sel);
}

void System_atoms_dup2(System* s, boost::python::list& data){
    System_atoms_dup1(s,data,NULL);
}

void System_atoms_add1(System* s, boost::python::list& atm,
                       boost::python::list& crd,
                       Selection* sel){
    vector<Atom> a;
    vector<Vector3f> c;
    a.resize(len(atm));
    c.resize(len(crd));
    for(int i=0;i<a.size();++i) a[i] = extract<Atom>(atm[i]);
    for(int i=0;i<c.size();++i){
        boost::python::object o = crd[i];
        MAP_EIGEN_TO_PYARRAY(v,Vector3f,o.ptr())
        c[i] = v;
    }
    s->atoms_add(a,c,sel);
}

void System_atoms_add2(System* s, boost::python::list& atm, boost::python::list& crd){
    System_atoms_add1(s,atm,crd,NULL);
}

float System_distance1(System* s, int i, int j, int fr, bool is_periodic){
    s->distance(i,j,fr,is_periodic);
}

float System_distance2(System* s, int i, int j, int fr){
    s->distance(i,j,fr,false);
}

float System_distance3(System* s, PyObject* p1, PyObject* p2, int fr, bool is_periodic){
    if(PyArray_Check(p1)){
        MAP_EIGEN_TO_PYARRAY(_p1,Vector3f,p1)
        MAP_EIGEN_TO_PYARRAY(_p2,Vector3f,p2)
        s->distance(_p1,_p2,fr,is_periodic);
    } else {
        System_distance1(s,extract<int>(p1),extract<int>(p2),fr,is_periodic);
    }
}

float System_distance4(System* s, PyObject* p1, PyObject* p2, int fr){
    System_distance3(s,p1,p2,fr,false);
}

void System_wrap_to_box(System* s, int frame, PyObject* point){
    MAP_EIGEN_TO_PYARRAY(p,Vector3f,point)
    Vector3f pp;
    s->wrap_to_box(frame,pp);
    p = pp;
}

PyObject* System_get_closest_image(System* s, PyObject* point, PyObject* target, int fr){
    MAP_EIGEN_TO_PYARRAY(p,Vector3f,point)
    MAP_EIGEN_TO_PYARRAY(t,Vector3f,target)
    Vector3f pp,tt;
    pp = p;
    tt = t;
    CREATE_PYARRAY_1D(ret,3)
    MAP_EIGEN_TO_PYARRAY(r,Vector3f,ret)
    r = s->get_closest_image(pp,tt,fr);
    return boost::python::incref(ret);
}

//-------------------------
// For Frame
//-------------------------

boost::python::list Frame_get_coord(Frame* f){
    boost::python::list l;
    for(int i=0;i<f->coord.size();++i){
        CREATE_PYARRAY_1D(p,3)
        MAP_EIGEN_TO_PYARRAY(v,Vector3f,p)
        v = f->coord[i];
        l.append(handle<>(p));
    }
    return l;
}

void Frame_set_coord(Frame* f, boost::python::list l){
    f->coord.resize(len(l));
    for(int i=0;i<f->coord.size();++i){
        boost::python::object o = l[i];
        MAP_EIGEN_TO_PYARRAY(v,Vector3f,o.ptr())
        f->coord[i] = v;
    }
}

PyObject* Frame_get_box(Frame* f){
    CREATE_PYARRAY_2D(p,3,3)
    MAP_EIGEN_TO_PYARRAY(m,Matrix3f,p)
    m = f->box;
    return boost::python::incref(p);
}
void Frame_set_box(Frame* f, PyObject* arr){
    MAP_EIGEN_TO_PYARRAY(m,Matrix3f,arr)
     f->box = m;
}


void make_bindings_System(){
    ///////////////////////////////////////////
    // Bindings for System
    ///////////////////////////////////////////
    import_array();

    class_<Frame>("Frame", init<>())
        .add_property("coord",&Frame_get_coord,&Frame_set_coord)
        .def_readwrite("t", &Frame::t)
        .add_property("box",&Frame_get_box,&Frame_set_box)
    ;

    class_<System, boost::noncopyable>("System", init<>())
        .def(init<std::string>() )
        .def(init<System>() )
        .def("num_atoms", &System::num_atoms)
        .def("num_frames", &System::num_frames)
        .def("load", &System::load, load_overloads())
        .def("frame_dup", &System::frame_dup)
        .def("set_frame", &System::set_frame)
        .def("frame_copy", &System::frame_copy)
        .def("frame_delete", &System::frame_delete, frame_delete_overloads())
        .def("getFrame_data", &System_getFrame_data)
        .def("setFrame_data", &System_setFrame_data)
        .def("update_selections", &System::update_selections)
        .def("getBox", &System_getBox)
        .def("setBox", &System_setBox)
        .def("is_box_triclinic", &System::is_box_triclinic)
        .def("getTime", &System_getTime)
        .def("setTime", &System_setTime)
        .def("get_box_vectors_angles",&System_get_box_vectors_angles)        
        .def("getXYZ", &System_getXYZ)
        .def("setXYZ", &System_setXYZ)
        .def("frame_append", &System::frame_append)
        .def("assign_resindex", &System::assign_resindex)
        .def("atoms_dup", &System_atoms_dup1)
        .def("atoms_dup", &System_atoms_dup2)
        .def("atoms_add", &System_atoms_add1)
        .def("atoms_add", &System_atoms_add2)
        .def("distance", &System_distance1)
        .def("distance", &System_distance2)
        .def("distance", &System_distance3)
        .def("distance", &System_distance4)
        .def("wrap_to_box", &System_wrap_to_box)
        .def("get_closest_image", &System_get_closest_image)
    ;
}
