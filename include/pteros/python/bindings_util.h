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

#ifndef BINDINGS_UTIL_H
#define BINDINGS_UTIL_H

#include "pteros/core/pteros_error.h"
#include "pteros/analysis/trajectory_processor.h"

#include <boost/python.hpp>
#include <numpy/noprefix.h>

#include <boost/python/numeric.hpp>

#include <iostream>

namespace pteros {

// Some macros

#define MAP_EIGEN_TO_PYARRAY(T, matr, pyobj) \
    if(!PyArray_Check(pyobj)) throw pteros::Pteros_error("NumPy array expected!"); \
    if(PyArray_TYPE(pyobj)!=PyArray_FLOAT) throw pteros::Pteros_error("float NumPy array expected!"); \
    Eigen::Map<T> matr((float*) PyArray_DATA(pyobj), \
            (PyArray_DIM((PyArrayObject*)pyobj,0)==PyArray_Size(pyobj)) ? PyArray_DIM((PyArrayObject*)pyobj,0) : PyArray_DIM((PyArrayObject*)pyobj,1), \
            (PyArray_DIM((PyArrayObject*)pyobj,0)==PyArray_Size(pyobj)) ? 1 : PyArray_DIM((PyArrayObject*)pyobj,0) );


template<class T>
void mapper(PyObject* pyobj, Eigen::Map<T>* map_ptr, PyObject* aux){
    size_t py_dim1, py_dim2, dim1, dim2;
    // Check if pyobj is a numpy array
    if(PyArray_Check(pyobj)){
        py_dim1 = PyArray_DIM((PyArrayObject*)pyobj,0);
        py_dim2 = PyArray_DIM((PyArrayObject*)pyobj,1);
        dim1 = (PyArray_Size(pyobj)==py_dim1) ? py_dim1 : py_dim2;
        dim2 = (PyArray_Size(pyobj)==py_dim1) ? 1 : py_dim1;
        // For fixed size Eigen object do a sanity check
        if(    T::RowsAtCompileTime!=Eigen::Dynamic
            && T::ColsAtCompileTime!=Eigen::Dynamic
            && (T::RowsAtCompileTime!=dim1 || T::ColsAtCompileTime!=dim2) )
        {
            throw Pteros_error()<<"Passed an object of size "<<dim1<<":"<<dim2
                            << " while expected size is "
                            << T::RowsAtCompileTime<<":"<<T::ColsAtCompileTime<<"!";
        }

        // Check array type
        if(PyArray_TYPE(pyobj)==PyArray_FLOAT){
            // The best case! Direct mapping
            new (map_ptr) Eigen::Map<T>((float*) PyArray_DATA(pyobj),dim1,dim2);
        } else {
            // Wrong type :( Need to cast.
            aux = PyArray_Cast((PyArrayObject*)pyobj,PyArray_FLOAT);
            new (map_ptr) Eigen::Map<T>((float*) PyArray_DATA(aux),dim1,dim2);
        }

    } else if(PySequence_Check(pyobj)) { // Check if it is a sequence
        // Get first dimension as list size
        py_dim1 = PySequence_Size(pyobj);
        // Get second dimension as the size of first list element
        PyObject* item0 = PySequence_GetItem(pyobj,0);
        if(PySequence_Check(item0)){
            py_dim2 = PySequence_Size(item0);
        } else {
            py_dim2 = 1;
        }

        dim1 = (py_dim2==1) ? py_dim1 : py_dim2;
        dim2 = (py_dim2==1) ? 1 : py_dim1;

        // For fixed size Eigen object do a sanity check
        if(    T::RowsAtCompileTime!=Eigen::Dynamic
            && T::ColsAtCompileTime!=Eigen::Dynamic
            && (T::RowsAtCompileTime!=dim1 || T::ColsAtCompileTime!=dim2) )
        {
            throw Pteros_error()<< "Passed an object of size "<<dim1<<":"<<dim2
                            << " while expected size is "
                            << T::RowsAtCompileTime<<":"<<T::ColsAtCompileTime<<"!";
        }
        // Convert to array and map
        aux = PyArray_FROM_OT(pyobj,PyArray_FLOAT);
        new (map_ptr) Eigen::Map<T>((float*) PyArray_DATA(aux),dim1,dim2);
    }
}

#define MAP_EIGEN_TO_PYTHON(T,matr,pyobj) \
    Eigen::Map<T> matr(NULL); \
    PyObject* _aux_for_##matr; \
    mapper(pyobj,&matr,_aux_for_##matr);

#define CREATE_PYARRAY_1D(pyobj, dim) \
    PyObject* pyobj; \
    { \
        npy_intp sz[1] = {dim};\
        pyobj = PyArray_SimpleNew(1, sz, PyArray_FLOAT); \
    }

#define CREATE_PYARRAY_1D_AND_MAP(pyobj, T, matr, dim) \
    PyObject* pyobj; \
    { \
        npy_intp sz[1] = {dim};\
        pyobj = PyArray_SimpleNew(1, sz, PyArray_FLOAT); \
    } \
    Eigen::Map<T> matr((float*)PyArray_DATA(pyobj),1,dim);


#define CREATE_PYARRAY_2D(pyobj, dim1, dim2) \
    PyObject* pyobj; \
    { \
        npy_intp dims[2] = {dim2,dim1}; \
        pyobj = PyArray_SimpleNew(2, dims, PyArray_FLOAT); \
    }

#define CREATE_PYARRAY_2D_AND_MAP(pyobj, T, matr, dim1, dim2) \
    PyObject* pyobj; \
    { \
        npy_intp dims[2] = {dim2,dim1}; \
        pyobj = PyArray_SimpleNew(2, dims, PyArray_FLOAT); \
    } \
    Eigen::Map<T> matr((float*)PyArray_DATA(pyobj),dim1,dim2);


} // End of namespace Pteros

#endif
