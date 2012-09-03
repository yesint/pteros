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

#include "bindings_util.h"

namespace pteros {

// Aux functions to determine array dimensions
/*
int get_dim1_from_PyArray(PyObject* arr){
    int dim1 = PyArray_DIM((PyArrayObject*)arr,0);
    int dim2 = PyArray_DIM((PyArrayObject*)arr,1);
    if(dim1==PyArray_Size(arr)){
        // 1D array
        //R = _dim1;
        //C = 1;
        return dim1;
    } else {
        // 2D array
        //R = _dim2;
        //C = _dim1;
        return dim2;
    }
}

int get_dim2_from_PyArray(PyObject* arr){
    int dim1 = PyArray_DIM((PyArrayObject*)arr,0);
    //int dim2 = PyArray_DIM((PyArrayObject*)arr,1);
    if(dim1==PyArray_Size(arr)){
        // 1D array
        //R = _dim1;
        //C = 1;
        return 1;
    } else {
        // 2D array
        //R = _dim2;
        //C = _dim1;
        return dim1;
    }
}
*/
} // End of namespace Pteros

