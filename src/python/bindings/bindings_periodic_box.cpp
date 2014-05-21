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

#include "bindings_periodic_box.h"
#include "pteros/core/system.h"
#include "pteros/python/bindings_util.h"

using namespace pteros;
using namespace Eigen;
using namespace boost::python;

void Periodic_box_modify(Periodic_box* b, PyObject* arr){
    MAP_EIGEN_TO_PYTHON_F(Matrix3f,m,arr)
    b->modify(m);
}

PyObject* Periodic_box_get_box(Periodic_box* b){
    CREATE_PYARRAY_2D_AND_MAP(p,Matrix3f,m,3,3)
    m = b->get_box();
    return boost::python::incref(p);
}

PyObject* Periodic_box_to_box(Periodic_box* b, PyObject* point){
    CREATE_PYARRAY_1D_AND_MAP(ret,Vector3f,v,3)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p,point)
    v = b->lab_to_box(p);
    return boost::python::incref(ret);
}

PyObject* Periodic_box_to_lab(Periodic_box* b, PyObject* point){
    CREATE_PYARRAY_1D_AND_MAP(ret,Vector3f,v,3)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p,point)
    v = b->box_to_lab(p);
    return boost::python::incref(ret);
}

PyObject* Periodic_box_to_box_matrix(Periodic_box* b){
    CREATE_PYARRAY_2D_AND_MAP(ret,Matrix3f,m,3,3)
    m = b->lab_to_box_matrix();
    return boost::python::incref(ret);
}

PyObject* Periodic_box_to_lab_matrix(Periodic_box* b){
    CREATE_PYARRAY_2D_AND_MAP(ret,Matrix3f,m,3,3)
    m = b->box_to_lab_matrix();
    return boost::python::incref(ret);
}

PyObject* Periodic_box_extents(Periodic_box* b){
    CREATE_PYARRAY_1D_AND_MAP(ret,Vector3f,v,3)
    v = b->extents();
    return boost::python::incref(ret);
}

float Periodic_box_distance(Periodic_box* b, PyObject* point1, PyObject* point2,
                             bool do_wrap=true, PyObject* periodic_dims=nullptr){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p1,point1)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p2,point2)
    if(periodic_dims){
        MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,periodic_dims)
        return b->distance(p1,p2,do_wrap,dim);
    } else {
        return b->distance(p1,p2,do_wrap);
    }
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Periodic_box_distance_overloads,Periodic_box_distance,3,5);

void Periodic_box_wrap_point(Periodic_box* b, PyObject* point, PyObject* dims_to_wrap=nullptr){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p,point)
    if(dims_to_wrap){
        MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims_to_wrap)
        b->wrap_point(p,dim);
    } else {
        b->wrap_point(p);
    }
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Periodic_box_wrap_point_overloads,Periodic_box_wrap_point,2,3);

PyObject* Periodic_box_get_closest_image(Periodic_box* b, PyObject* point, PyObject* target,
                             bool do_wrap=true, PyObject* dims_to_wrap=nullptr){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p,point)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,t,target)   
    CREATE_PYARRAY_1D_AND_MAP(ret,Vector3f,v,3)
    if(dims_to_wrap){
        MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims_to_wrap)
        v = b->get_closest_image(p,t,do_wrap,dim);
        return boost::python::incref(ret);
    } else {
        v = b->get_closest_image(p,t,do_wrap);
        return boost::python::incref(ret);
    }
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Periodic_box_get_closest_image_overloads,Periodic_box_get_closest_image,3,5);

PyObject* Periodic_box_shortest_vector(Periodic_box* b, PyObject* point1, PyObject* point2,
                                        PyObject* dims_to_wrap=nullptr){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p1,point1)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p2,point2)    
    CREATE_PYARRAY_1D_AND_MAP(ret,Vector3f,v,3)
    if(dims_to_wrap){
        MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims_to_wrap)
        v = b->shortest_vector(p1,p2,dim);
        return boost::python::incref(ret);
    } else {
        v = b->shortest_vector(p1,p2);
        return boost::python::incref(ret);
    }
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Periodic_box_shortest_vector_overloads,Periodic_box_shortest_vector,3,4);


void make_bindings_Periodic_box(){
    import_array();

    class_<Periodic_box>("Periodic_box", init<>())
        .def("modify",&Periodic_box_modify)
        .def("get_box",&Periodic_box_get_box)
        .def("lab_to_box",&Periodic_box_to_box)
        .def("lab_to_box_matrix",&Periodic_box_to_box_matrix)
        .def("box_to_lab",&Periodic_box_to_lab)
        .def("box_to_lab_matrix",&Periodic_box_to_lab_matrix)
        .def("extent",&Periodic_box::extent)
        .def("extents",&Periodic_box_extents)
        .def("is_triclinic",&Periodic_box::is_triclinic)
        .def("is_periodic",&Periodic_box::is_periodic)
        .def("distance",&Periodic_box_distance, Periodic_box_distance_overloads())
        .def("wrap_point",&Periodic_box_wrap_point, Periodic_box_wrap_point_overloads())
        .def("get_closest_image",&Periodic_box_get_closest_image, Periodic_box_get_closest_image_overloads())
        .def("shortest_vector",&Periodic_box_shortest_vector, Periodic_box_shortest_vector_overloads())
        .def("is_periodic",&Periodic_box::volume)
    ;
}
