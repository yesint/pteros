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


#define CREATE_PYARRAY_1D(pyobj, dim) \
    PyObject* pyobj; \
    { \
        npy_intp sz[1] = {dim};\
        pyobj = PyArray_SimpleNew(1, sz, PyArray_FLOAT); \
    }

#define CREATE_PYARRAY_2D(pyobj, dim1, dim2) \
    PyObject* pyobj; \
    { \
        npy_intp dims[2] = {dim2,dim1}; \
        pyobj = PyArray_SimpleNew(2, dims, PyArray_FLOAT); \
    }

} // End of namespace Pteros

#endif
