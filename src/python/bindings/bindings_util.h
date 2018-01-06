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

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>


namespace py = pybind11;

// Aux function to create py::array on top of std::vector<POD>* which owns the vector
template<class T>
py::array vector_to_array(std::vector<T>* ptr, size_t sz=-1){
    if(sz==-1) sz=ptr->size();
    auto capsule = py::capsule(ptr, [](void *v) { delete reinterpret_cast<std::vector<T>*>(v); });
    return py::array(sz, ptr->data(), capsule);
}