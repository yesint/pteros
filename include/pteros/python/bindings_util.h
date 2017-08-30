/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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
#include <boost/python.hpp>
#include <numpy/noprefix.h>
#include <boost/python/numeric.hpp>
#include <iostream>
#include <Eigen/Core>

namespace pteros {

// Some helper functions and macros

// Class for releasing resources
class Decrementer {
public:
    Decrementer(PyObject* pyobj): m_obj(pyobj) {}
    ~Decrementer(){ Py_XDECREF(m_obj); }
private:
    PyObject* m_obj;
};

template<class T, int PyT>
PyObject* mapper(PyObject* pyobj, size_t& dim1, size_t& dim2){
    PyObject* aux = nullptr;    

    size_t py_dim1, py_dim2;

    // Check if pyobj is a numpy array
    if(PyArray_Check(pyobj)){                
        int ndim = PyArray_NDIM((PyArrayObject*)pyobj);
        if(ndim>2) throw Pteros_error("Need 1D or 2D array!");

        py_dim1 = PyArray_DIM((PyArrayObject*)pyobj,0);
        py_dim2 = PyArray_DIM((PyArrayObject*)pyobj,1);
        dim1 = (ndim==1) ? py_dim1 : py_dim2;
        dim2 = (ndim==1) ? 1 : py_dim1;

        // For fixed size Eigen object do a sanity check
        if(    T::RowsAtCompileTime!=Eigen::Dynamic
            && T::ColsAtCompileTime!=Eigen::Dynamic
            && (T::RowsAtCompileTime!=dim1 || T::ColsAtCompileTime!=dim2) )
        {
            throw Pteros_error("Passed an array of size {}:{} while expected size is {}:{}",
                               dim1, dim2, T::RowsAtCompileTime, T::ColsAtCompileTime);
        }

        // Check array type
        if(PyArray_TYPE((PyArrayObject*)pyobj)==PyT){
            // The best case! Direct mapping
            // We need to increase reference counter of initial object
            // since it will be decreased later by Decrementer!
            Py_XINCREF((PyObject*)pyobj);
            return (PyObject*)pyobj;
        } else {
            // Wrong type :( Need to cast.
            aux = PyArray_Cast((PyArrayObject*)pyobj,PyT);
            return aux;
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
        Py_XDECREF(item0);

        dim1 = (py_dim2==1) ? py_dim1 : py_dim2;
        dim2 = (py_dim2==1) ? 1 : py_dim1;

        // For fixed size Eigen object do a sanity check
        if(    T::RowsAtCompileTime!=Eigen::Dynamic
            && T::ColsAtCompileTime!=Eigen::Dynamic
            && (T::RowsAtCompileTime!=dim1 || T::ColsAtCompileTime!=dim2) )
        {
            throw Pteros_error("Passed an object of size {}:{} while expected size is {}:{}",
                               dim1, dim2, T::RowsAtCompileTime, T::ColsAtCompileTime);
        }

        // Convert to array and map
        aux = PyArray_FROM_OT(pyobj,PyT);
        if(!aux) throw Pteros_error("Can't convert sequence to array!");
        return aux;
    } else {
        // Something incompatible
        throw Pteros_error("An array or a sequence is required!");
    }
}

#define MAP_EIGEN_TO_PYTHON_F(T,matr,pyobj) \
    size_t __dim1__for__##matr, __dim2__for__##matr; \
    PyObject* __pyobj__for__##matr = mapper<T,NPY_FLOAT>(pyobj,__dim1__for__##matr,__dim2__for__##matr);\
    Decrementer __dec__for__##matr(__pyobj__for__##matr);\
    Eigen::Map<T>  matr((float*)PyArray_DATA((PyArrayObject*)__pyobj__for__##matr),__dim1__for__##matr,__dim2__for__##matr);

#define MAP_EIGEN_TO_PYTHON_I(T,matr,pyobj) \
    size_t __dim1__for__##matr, __dim2__for__##matr; \
    PyObject* __pyobj__for__##matr = mapper<T,NPY_INT>(pyobj,__dim1__for__##matr,__dim2__for__##matr);\
    Decrementer __dec__for__##matr(__pyobj__for__##matr);\
    Eigen::Map<T>  matr((int*)PyArray_DATA((PyArrayObject*)__pyobj__for__##matr),__dim1__for__##matr,__dim2__for__##matr);


#define CREATE_PYARRAY_1D_AND_MAP_F(pyobj, T, matr, dim) \
    PyObject* pyobj; \
    { \
        npy_intp sz[1] = {dim};\
        pyobj = PyArray_SimpleNew(1, sz, NPY_FLOAT); \
    } \
    Eigen::Map<T> matr((float*)PyArray_DATA((PyArrayObject*)pyobj),dim,1);

#define CREATE_PYARRAY_1D_AND_MAP_I(pyobj, T, matr, dim) \
    PyObject* pyobj; \
    { \
        npy_intp sz[1] = {dim};\
        pyobj = PyArray_SimpleNew(1, sz, NPY_INT); \
    } \
    Eigen::Map<T> matr((int*)PyArray_DATA((PyArrayObject*)pyobj),dim,1);


#define CREATE_PYARRAY_2D_AND_MAP_F(pyobj, T, matr, dim1, dim2) \
    PyObject* pyobj; \
    { \
        npy_intp dims[2] = {dim2,dim1}; \
        pyobj = PyArray_SimpleNew(2, dims, NPY_FLOAT); \
    } \
    Eigen::Map<T> matr((float*)PyArray_DATA((PyArrayObject*)pyobj),dim1,dim2);

#define CREATE_PYARRAY_2D_AND_MAP_I(pyobj, T, matr, dim1, dim2) \
    PyObject* pyobj; \
    { \
        npy_intp dims[2] = {dim2,dim1}; \
        pyobj = PyArray_SimpleNew(2, dims, NPY_INT); \
    } \
    Eigen::Map<T> matr((int*)PyArray_DATA((PyArrayObject*)pyobj),dim1,dim2);

#define PYOBJECT_TO_VECTOR(pyobj,T,vecname) \
    vector<T> vecname; \
    if( PySequence_Check(pyobj) ){ \
        boost::python::object __list__for__##pyobj( handle<>(borrowed(pyobj)) ); \
        int n = boost::python::len(__list__for__##pyobj); \
        vecname.resize(n); \
        for(int i=0;i<n;++i) vecname[i] = extract<T>(__list__for__##pyobj[i]); \
    } else { \
        throw Pteros_error("Provided python object is not a sequence!"); \
    }


} // End of namespace Pteros

#endif
