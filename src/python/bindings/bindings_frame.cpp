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

struct Frame_suite : boost::python::pickle_suite
{
    static boost::python::tuple getstate(boost::python::object w_obj)
    {
        Frame const& w = boost::python::extract<Frame const&>(w_obj)();
        return boost::python::make_tuple(
                  w_obj.attr("__dict__"), //If the python object has other attributes, they will be stored in the dict
                  w.box,
                  w.coord,
                  w.time);
    }

    static void setstate(boost::python::object w_obj, boost::python::tuple state)
    {
        using namespace boost::python;
        Frame& w = extract<Frame&>(w_obj)();
        // restore the object's __dict__
        dict d = extract<dict>(w_obj.attr("__dict__"))();
        d.update(state[0]);
        //w.box = extract<int>(state[1]);
        //w.coord = extract<float>(state[2]);
        //w.time = extract<int>(state[3]);
    }
    static bool getstate_manages_dict() { return true; }
};

boost::python::list Frame_get_coord(Frame* f){
    boost::python::list l;
    for(int i=0;i<f->coord.size();++i){
        CREATE_PYARRAY_1D_AND_MAP_F(p,Vector3f,v,3)
        v = f->coord[i];
        l.append(handle<>(p));
    }
    return l;
}

PyObject* Frame_get_coord_array(Frame* f){
    CREATE_PYARRAY_2D_AND_MAP_F(p,MatrixXf,m,3,npy_intp(f->coord.size()))
    for(int i=0;i<f->coord.size();++i){
        m.col(i) = f->coord[i];
    }
    return p;
}

void Frame_set_coord_array(Frame* f, PyObject* arr){
    MAP_EIGEN_TO_PYTHON_F(MatrixXf,m,arr)
    f->coord.reserve(m.cols());
    for(int i=0;i<m.cols();++i){
        f->coord.push_back(m.col(i));
    }
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
        .enable_pickling()
        .def_pickle(Frame_suite())
        .def("get_coord_array",&Frame_get_coord_array)
        .def("set_coord_array",&Frame_set_coord_array)
    ;
}
