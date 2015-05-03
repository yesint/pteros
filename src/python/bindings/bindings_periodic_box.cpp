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

boost::shared_ptr<Periodic_box> constructor0(PyObject* vectors, PyObject* angles){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,v,vectors)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,a,angles)
    boost::shared_ptr<Periodic_box> g(new Periodic_box(v,a));
    return g;
}

void Periodic_box_modify(Periodic_box* b, PyObject* arr){
    MAP_EIGEN_TO_PYTHON_F(Matrix3f,m,arr)
    b->modify(m);
}

PyObject* Periodic_box_get_matrix(Periodic_box* b){
    CREATE_PYARRAY_2D_AND_MAP_F(p,Matrix3f,m,3,3)
    m = b->get_matrix();
    return boost::python::incref(p);
}

PyObject* Periodic_box_get_inv_matrix(Periodic_box* b){
    CREATE_PYARRAY_2D_AND_MAP_F(p,Matrix3f,m,3,3)
    m = b->get_inv_matrix();
    return boost::python::incref(p);
}

PyObject* Periodic_box_get_vector(Periodic_box* b, int i){
    CREATE_PYARRAY_1D_AND_MAP_F(vec,Vector3f,v,3)
    v = b->get_vector(i);
    return boost::python::incref(vec);
}

void Periodic_box_scale_vectors(Periodic_box* b, PyObject* scale){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,s,scale)
    b->scale_vectors(s);
}

PyObject* Periodic_lab_to_box(Periodic_box* b, PyObject* point){
    CREATE_PYARRAY_1D_AND_MAP_F(ret,Vector3f,v,3)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p,point)
    v = b->lab_to_box(p);
    return boost::python::incref(ret);
}

PyObject* Periodic_box_to_lab(Periodic_box* b, PyObject* point){
    CREATE_PYARRAY_1D_AND_MAP_F(ret,Vector3f,v,3)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p,point)
    v = b->box_to_lab(p);
    return boost::python::incref(ret);
}

PyObject* Periodic_lab_to_box_transform(Periodic_box* b){
    CREATE_PYARRAY_2D_AND_MAP_F(ret,Matrix3f,m,3,3)
    m = b->lab_to_box_transform();
    return boost::python::incref(ret);
}

PyObject* Periodic_box_to_lab_transform(Periodic_box* b){
    CREATE_PYARRAY_2D_AND_MAP_F(ret,Matrix3f,m,3,3)
    m = b->box_to_lab_transform();
    return boost::python::incref(ret);
}

PyObject* Periodic_box_extents(Periodic_box* b){
    CREATE_PYARRAY_1D_AND_MAP_F(ret,Vector3f,v,3)
    v = b->extents();
    return boost::python::incref(ret);
}

float Periodic_box_distance(Periodic_box* b, PyObject* point1, PyObject* point2,
                             PyObject* periodic_dims=nullptr){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p1,point1)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p2,point2)
    if(periodic_dims){
        MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,periodic_dims)
        return b->distance(p1,p2,dim);
    } else {
        return b->distance(p1,p2);
    }
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Periodic_box_distance_overloads,Periodic_box_distance,3,4);

float Periodic_box_distance_squared(Periodic_box* b, PyObject* point1, PyObject* point2,
                             PyObject* periodic_dims=nullptr){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p1,point1)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p2,point2)
    if(periodic_dims){
        MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,periodic_dims)
        return b->distance_squared(p1,p2,dim);
    } else {
        return b->distance_squared(p1,p2);
    }
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Periodic_box_distance_squared_overloads,Periodic_box_distance_squared,3,4);

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
                             PyObject* dims_to_wrap=nullptr){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p,point)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,t,target)   
    CREATE_PYARRAY_1D_AND_MAP_F(ret,Vector3f,v,3)
    if(dims_to_wrap){
        MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims_to_wrap)
        v = b->get_closest_image(p,t,dim);
        return boost::python::incref(ret);
    } else {
        v = b->get_closest_image(p,t);
        return boost::python::incref(ret);
    }
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Periodic_box_get_closest_image_overloads,Periodic_box_get_closest_image,3,4);

PyObject* Periodic_box_shortest_vector(Periodic_box* b, PyObject* point1, PyObject* point2,
                                        PyObject* dims_to_wrap=nullptr){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p1,point1)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p2,point2)    
    CREATE_PYARRAY_1D_AND_MAP_F(ret,Vector3f,v,3)
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


bool Periodic_box_in_box(Periodic_box* b, PyObject* point){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,p,point)
    return b->in_box(p);
}


boost::python::tuple Periodic_box_to_vectors_angles(Periodic_box* b){
    CREATE_PYARRAY_1D_AND_MAP_F(vec,Vector3f,v,3)
    CREATE_PYARRAY_1D_AND_MAP_F(ang,Vector3f,a,3)
    b->to_vectors_angles(v,a);
    return boost::python::make_tuple(handle<>(vec),handle<>(ang));
}

void Periodic_box_from_vectors_angles(Periodic_box* b, PyObject* vec, PyObject* ang){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,v,vec)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,a,ang)
    b->from_vectors_angles(v,a);
}


void make_bindings_Periodic_box(){
    import_array();

    class_<Periodic_box>("Periodic_box", init<>())
        .def("__init__",make_constructor(&constructor0))
        .def("modify",&Periodic_box_modify)
        .def("get_matrix",&Periodic_box_get_matrix)
        .def("get_inv_matrix",&Periodic_box_get_inv_matrix)
        .def("get_vector",&Periodic_box_get_vector)
        .def("scale_vectors",&Periodic_box_scale_vectors)
        .def("lab_to_box",&Periodic_lab_to_box)
        .def("lab_to_box_transform",&Periodic_lab_to_box_transform)
        .def("box_to_lab",&Periodic_box_to_lab)
        .def("box_to_lab_transform",&Periodic_box_to_lab_transform)
        .def("extent",&Periodic_box::extent)
        .def("extents",&Periodic_box_extents)
        .def("is_triclinic",&Periodic_box::is_triclinic)
        .def("is_periodic",&Periodic_box::is_periodic)
        .def("distance",&Periodic_box_distance, Periodic_box_distance_overloads())
        .def("distance_squared",&Periodic_box_distance_squared, Periodic_box_distance_squared_overloads())
        .def("wrap_point",&Periodic_box_wrap_point, Periodic_box_wrap_point_overloads())
        .def("in_box",&Periodic_box_in_box)
        .def("volume",&Periodic_box::volume)
        .def("get_closest_image",&Periodic_box_get_closest_image, Periodic_box_get_closest_image_overloads())
        .def("shortest_vector",&Periodic_box_shortest_vector, Periodic_box_shortest_vector_overloads())
        .def("to_vectors_angles",&Periodic_box_to_vectors_angles)
        .def("from_vectors_angles",&Periodic_box_from_vectors_angles)
    ;
}
