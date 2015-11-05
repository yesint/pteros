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

#include "bindings_frame.h"
#include "pteros/core/system.h"
#include "pteros/python/bindings_util.h"

using namespace pteros;
using namespace Eigen;
using namespace boost::python;

boost::python::list Frame_get_coord(Frame* f){
    boost::python::list l;
    for(int i=0;i<f->coord.size();++i){
        CREATE_PYARRAY_1D_AND_MAP_F(p,Vector3f,v,3)
        v = f->coord[i];
        l.append(handle<>(p));
    }
    return l;
}

void Frame_set_coord(Frame* f, boost::python::list l){
    f->coord.resize(len(l));
    for(int i=0;i<f->coord.size();++i){
        boost::python::object o = l[i];
        MAP_EIGEN_TO_PYTHON_F(Vector3f,v,o.ptr())
        f->coord[i] = v;
    }
}

PyObject* Frame_get_box(Frame* f){
    CREATE_PYARRAY_2D_AND_MAP_F(p,Matrix3f,m,3,3)
    m = f->box.get_matrix();
    return p;
}

void Frame_set_box(Frame* f, PyObject* arr){
    MAP_EIGEN_TO_PYTHON_F(Matrix3f,m,arr)
    f->box.modify(m);
}

void make_bindings_Frame(){
    import_array();

    class_<Frame>("Frame", init<>())
        .add_property("coord",&Frame_get_coord,&Frame_set_coord)
        .def_readwrite("t", &Frame::time)
        .add_property("box",&Frame_get_box,&Frame_set_box)
    ;
}
