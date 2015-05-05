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

#include "bindings_utilities.h"
#include "pteros/python/bindings_util.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>

using namespace boost::python;
using namespace pteros;
using namespace Eigen;

float angle_between_vectors_py(PyObject* vec1, PyObject* vec2){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,v1,vec1)
    MAP_EIGEN_TO_PYTHON_F(Vector3f,v2,vec2)
    return angle_between_vectors(v1,v2);
}

void make_bindings_utilities(){
    import_array();
    def("angle_between_vectors",&angle_between_vectors_py);
}
