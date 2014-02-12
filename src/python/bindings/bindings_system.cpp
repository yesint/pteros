/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2014, Semen Yesylevskyy
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

using namespace pteros;
using namespace Eigen;
using namespace boost::python;

/**********************
  Wrappers for System
***********************/

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(load_overloads, load, 1, 3)
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

/*
float System_distance3(System* s, PyObject* p1, PyObject* p2, int fr, bool is_periodic){
    if(PyArray_Check(p1)){
        MAP_EIGEN_TO_PYARRAY(_p1,Vector3f,p1)
        MAP_EIGEN_TO_PYARRAY(_p2,Vector3f,p2)
        return s->Box(fr).distance(_p1,_p2);
    } else {
        return System_distance1(s,extract<int>(p1),extract<int>(p2),fr,is_periodic);
    }
}
*/

/*
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
*/

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
    m = f->box.get_box();
    return boost::python::incref(p);
}

void Frame_set_box(Frame* f, PyObject* arr){
    MAP_EIGEN_TO_PYARRAY(m,Matrix3f,arr)
     f->box.modify(m);
}

void Periodic_box_modify(Periodic_box* b, PyObject* arr){
    MAP_EIGEN_TO_PYARRAY(m,Matrix3f,arr)
    b->modify(m);
}

PyObject* Periodic_box_get_box(Periodic_box* b){
    CREATE_PYARRAY_2D(p,3,3)
    MAP_EIGEN_TO_PYARRAY(m,Matrix3f,p)
    m = b->get_box();
    return boost::python::incref(p);
}

PyObject* Periodic_box_to_box(Periodic_box* b, PyObject* point){
    CREATE_PYARRAY_1D(ret,3)
    MAP_EIGEN_TO_PYARRAY(v,Vector3f,ret)
    MAP_EIGEN_TO_PYARRAY(p,Vector3f,point)
    v = b->to_box(p);
    return boost::python::incref(ret);
}

PyObject* Periodic_box_to_lab(Periodic_box* b, PyObject* point){
    CREATE_PYARRAY_1D(ret,3)
    MAP_EIGEN_TO_PYARRAY(v,Vector3f,ret)
    MAP_EIGEN_TO_PYARRAY(p,Vector3f,point)
    v = b->to_lab(p);
    return boost::python::incref(ret);
}

PyObject* Periodic_box_to_box_matrix(Periodic_box* b){
    CREATE_PYARRAY_2D(ret,3,3)
    MAP_EIGEN_TO_PYARRAY(m,Matrix3f,ret)
    m = b->to_box_matrix();
    return boost::python::incref(ret);
}

PyObject* Periodic_box_to_lab_matrix(Periodic_box* b){
    CREATE_PYARRAY_2D(ret,3,3)
    MAP_EIGEN_TO_PYARRAY(m,Matrix3f,ret)
    m = b->to_lab_matrix();
    return boost::python::incref(ret);
}

PyObject* Periodic_box_extents(Periodic_box* b){
    CREATE_PYARRAY_1D(ret,3)
    MAP_EIGEN_TO_PYARRAY(v,Vector3f,ret)
    v = b->extents();
    return boost::python::incref(ret);
}

float Periodic_box_distance1(Periodic_box* b, PyObject* point1, PyObject* point2,
                             bool do_wrap, boost::python::list periodic_dims){
    MAP_EIGEN_TO_PYARRAY(p1,Vector3f,point1)
    MAP_EIGEN_TO_PYARRAY(p2,Vector3f,point2)
    Vector3i dims;
    for(int i=0;i<3;++i) dims(i) = extract<int>(periodic_dims[i]);

    float dist = b->distance(p1,p2,do_wrap,dims);
    return dist;
}

float Periodic_box_distance2(Periodic_box* b, PyObject* point1, PyObject* point2, bool do_wrap){
    MAP_EIGEN_TO_PYARRAY(p1,Vector3f,point1)
    MAP_EIGEN_TO_PYARRAY(p2,Vector3f,point2)

    float dist = b->distance(p1,p2,do_wrap);
    return dist;
}

float Periodic_box_distance3(Periodic_box* b, PyObject* point1, PyObject* point2){
    MAP_EIGEN_TO_PYARRAY(p1,Vector3f,point1)
    MAP_EIGEN_TO_PYARRAY(p2,Vector3f,point2)

    float dist = b->distance(p1,p2);
    return dist;
}

void Periodic_box_wrap_point1(Periodic_box* b, PyObject* point, boost::python::list dims_to_wrap){
    MAP_EIGEN_TO_PYARRAY(p,Vector3f,point)

    Vector3i dims;
    for(int i=0;i<3;++i) dims(i) = extract<int>(dims_to_wrap[i]);

    Vector3f pp = p;
    b->wrap_point(pp,dims);
    p = pp;
}

void Periodic_box_wrap_point2(Periodic_box* b, PyObject* point){
    MAP_EIGEN_TO_PYARRAY(p,Vector3f,point)
    Vector3f pp = p;
    b->wrap_point(pp);
    p = pp;
}

PyObject* Periodic_box_get_closest_image1(Periodic_box* b, PyObject* point, PyObject* target,
                             bool do_wrap, boost::python::list dims_to_wrap){
    MAP_EIGEN_TO_PYARRAY(p,Vector3f,point)
    MAP_EIGEN_TO_PYARRAY(t,Vector3f,target)

    CREATE_PYARRAY_1D(ret,3)
    MAP_EIGEN_TO_PYARRAY(v,Vector3f,ret)

    Vector3i dims;
    for(int i=0;i<3;++i) dims(i) = extract<int>(dims_to_wrap[i]);

    v = b->get_closest_image(p,t,do_wrap,dims);
    return boost::python::incref(ret);
}

PyObject* Periodic_box_get_closest_image2(Periodic_box* b, PyObject* point, PyObject* target,
                             bool do_wrap){
    MAP_EIGEN_TO_PYARRAY(p,Vector3f,point)
    MAP_EIGEN_TO_PYARRAY(t,Vector3f,target)

    CREATE_PYARRAY_1D(ret,3)
    MAP_EIGEN_TO_PYARRAY(v,Vector3f,ret)

    v = b->get_closest_image(p,t,do_wrap);
    return boost::python::incref(ret);
}

PyObject* Periodic_box_get_closest_image3(Periodic_box* b, PyObject* point, PyObject* target){
    MAP_EIGEN_TO_PYARRAY(p,Vector3f,point)
    MAP_EIGEN_TO_PYARRAY(t,Vector3f,target)

    CREATE_PYARRAY_1D(ret,3)
    MAP_EIGEN_TO_PYARRAY(v,Vector3f,ret)

    v = b->get_closest_image(p,t);
    return boost::python::incref(ret);
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

    class_<Periodic_box>("Periodic_box", init<>())
        .def("modify",&Periodic_box_modify)
        .def("get_box",&Periodic_box_get_box)
        .def("to_box",&Periodic_box_to_box)
        .def("to_box_matrix",&Periodic_box_to_box_matrix)
        .def("to_lab",&Periodic_box_to_lab)
        .def("to_lab_matrix",&Periodic_box_to_lab_matrix)
        .def("extent",&Periodic_box::extent)
        .def("extents",&Periodic_box_extents)
        .def("is_triclinic",&Periodic_box::is_triclinic)
        .def("is_periodic",&Periodic_box::is_periodic)
        .def("distance",&Periodic_box_distance1)
        .def("distance",&Periodic_box_distance2)
        .def("distance",&Periodic_box_distance3)
        .def("wrap_point",&Periodic_box_wrap_point1)
        .def("wrap_point",&Periodic_box_wrap_point2)
        .def("get_closest_image",&Periodic_box_get_closest_image1)
        .def("get_closest_image",&Periodic_box_get_closest_image2)
        .def("get_closest_image",&Periodic_box_get_closest_image3)
        .def("is_periodic",&Periodic_box::volume)
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
        .def("getBox", &System_getBox,return_value_policy<reference_existing_object>())
        .def("setBox", &System_setBox)        
        .def("getTime", &System_getTime)
        .def("setTime", &System_setTime)
        // Returns tuple of (vectors,angles)        
        .def("getXYZ", &System_getXYZ)
        .def("setXYZ", &System_setXYZ)
        .def("frame_append", &System::frame_append)
        .def("assign_resindex", &System::assign_resindex)
        .def("atoms_dup", &System_atoms_dup1)
        .def("atoms_dup", &System_atoms_dup2)
        .def("atoms_add", &System_atoms_add1)
        .def("atoms_add", &System_atoms_add2)             
        .def("wrap_all", &System::wrap_all)
    ;
}
