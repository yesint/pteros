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

#ifndef BINDINGS_UTIL_H
#define BINDINGS_UTIL_H

#include "pteros/core/pteros_error.h"
#include "pteros/analysis/trajectory_processor.h"

#include <boost/python.hpp>
#include <numpy/noprefix.h>

#include <boost/python/numeric.hpp>

#include <iostream>

using namespace std;
using namespace boost;
using namespace boost::python;
using namespace boost::python::numeric;
using namespace Eigen;
using namespace pteros;

namespace pteros {

// Some macros

#define MAP_EIGEN_TO_PYARRAY(_matr,_T,_obj_ptr) \
    if(!PyArray_Check(_obj_ptr)) throw Pteros_error("NumPy array expected!"); \
    if(PyArray_TYPE(_obj_ptr)!=PyArray_FLOAT) throw Pteros_error("float NumPy array expected!"); \
    Map<_T> _matr((float*) PyArray_DATA(_obj_ptr), \
            (PyArray_DIM((PyArrayObject*)_obj_ptr,0)==PyArray_Size(_obj_ptr)) ? PyArray_DIM((PyArrayObject*)_obj_ptr,0) : PyArray_DIM((PyArrayObject*)_obj_ptr,1), \
            (PyArray_DIM((PyArrayObject*)_obj_ptr,0)==PyArray_Size(_obj_ptr)) ? 1 : PyArray_DIM((PyArrayObject*)_obj_ptr,0) );


#define CREATE_PYARRAY_1D(_ptr_obj, _dim1) \
    PyObject* _ptr_obj; \
    { \
        npy_intp _sz_dim1[1];\
        _sz_dim1[0] = _dim1; \
        _ptr_obj = PyArray_SimpleNew(1, _sz_dim1, PyArray_FLOAT); \
    }

#define CREATE_PYARRAY_2D(_ptr_obj, _dim1, _dim2) \
    PyObject* _ptr_obj; \
    { \
        npy_intp _dims[2]; \
        _dims[1] = _dim1; \
        _dims[0] = _dim2; \
        _ptr_obj = PyArray_SimpleNew(2, _dims, PyArray_FLOAT); \
    }

#define CREATE_PYARRAY_2D_FROM_DATA(_ptr_obj, _dim1, _dim2, _data_ptr) \
    PyObject* _ptr_obj; \
    { \
        npy_intp _dims[2]; \
        _dims[1] = _dim1; \
        _dims[0] = _dim2; \
        _ptr_obj = PyArray_SimpleNewFromData(2, _dims, PyArray_FLOAT, (void*)_data_ptr); \
    }

} // End of namespace Pteros

#endif
