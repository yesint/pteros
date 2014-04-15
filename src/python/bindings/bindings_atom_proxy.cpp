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

#include "bindings_atom_proxy.h"
#include "pteros/core/selection.h"
#include "pteros/python/bindings_util.h"

using namespace pteros;
using namespace Eigen;
using namespace boost::python;

//-------------------------------------------------------

PyObject* Atom_proxy_getXYZ1(Atom_proxy* s){    
    CREATE_PYARRAY_1D(p,3)
    MAP_EIGEN_TO_PYARRAY(data,Vector3f,p)
    data = s->XYZ();
    return boost::python::incref(p);
}

PyObject* Atom_proxy_getXYZ2(Atom_proxy* s, int fr){
    CREATE_PYARRAY_1D(p,3)
    MAP_EIGEN_TO_PYARRAY(data,Vector3f,p)
    data = s->XYZ(fr);
    return boost::python::incref(p);
}

void Atom_proxy_setXYZ1(Atom_proxy* s, PyObject* obj){
    MAP_EIGEN_TO_PYARRAY(data,Vector3f,obj)
    s->XYZ() = data;
}

void Atom_proxy_setXYZ2(Atom_proxy* s, int fr, PyObject* obj){
    MAP_EIGEN_TO_PYARRAY(data,Vector3f,obj)
    s->XYZ(fr) = data;
}


// Macros to wrap an inline accessor function
#define WRAP_ACCESSOR_FR(_out_type, _func) \
    _out_type Atom_proxy_get##_func##1(Atom_proxy* s){ return s->_func(); } \
    void Atom_proxy_set##_func##1(Atom_proxy* s, _out_type val){ s->_func() = val; } \
    _out_type Atom_proxy_get##_func##2(Atom_proxy* s, int fr){ return s->_func(fr); } \
    void Atom_proxy_set##_func##2(Atom_proxy* s, int fr, _out_type val){ s->_func(fr) = val; }

#define WRAP_ACCESSOR(_out_type, _func) \
    _out_type Atom_proxy_get##_func(Atom_proxy* s){ return s->_func(); } \
    void Atom_proxy_set##_func(Atom_proxy* s, _out_type val){ s->_func() = val; }

#define DEF_PROPERTY(_prop,_func) \
    .add_property(#_prop, &Atom_proxy_get##_func, &Atom_proxy_set##_func)


WRAP_ACCESSOR_FR(float,X)
WRAP_ACCESSOR_FR(float,Y)
WRAP_ACCESSOR_FR(float,Z)
WRAP_ACCESSOR(int,Type)
WRAP_ACCESSOR(std::string,Type_name)
WRAP_ACCESSOR(std::string,Resname)
WRAP_ACCESSOR(char,Chain)
WRAP_ACCESSOR(std::string,Name)
WRAP_ACCESSOR(float,Mass)
WRAP_ACCESSOR(float,Charge)
WRAP_ACCESSOR(float,Beta)
WRAP_ACCESSOR(float,Occupancy)
WRAP_ACCESSOR(long,Resid)
WRAP_ACCESSOR(std::string,Tag)
WRAP_ACCESSOR(int,Resindex)
WRAP_ACCESSOR(int,Index)

void make_bindings_Atom_proxy(){
    import_array();

    class_<Atom_proxy>("Atom_proxy", init<>())
        .def(init<Selection*,int>() )

        // For coordinate accessors we should use setX instead of just X in Python
        // This is because Python don't respect void in return - all functions
        // with equal number of argumets are considered equivalent thus
        // "float X(int,int)" and "void X(int,float)" become the same function...
        // For non-coordinate accessors this is not needed but used to be consistent                
        .def("getX",&Atom_proxy_getX2)        
        .def("setX",&Atom_proxy_setX2)
        .add_property("x", &Atom_proxy_getX1, &Atom_proxy_setX1)

        .def("getY",&Atom_proxy_getY2)        
        .def("setY",&Atom_proxy_setY2)
        .add_property("y", &Atom_proxy_getY1, &Atom_proxy_setY1)

        .def("getZ",&Atom_proxy_getZ2)        
        .def("setZ",&Atom_proxy_setZ2)
        .add_property("z", &Atom_proxy_getZ1, &Atom_proxy_setZ1)

        .def("getXYZ",&Atom_proxy_getXYZ2)        
        .def("setXYZ",&Atom_proxy_setXYZ2)
        .add_property("xyz", &Atom_proxy_getXYZ1, &Atom_proxy_setXYZ1)

        DEF_PROPERTY(type,Type)
        DEF_PROPERTY(type_name,Type_name)
        DEF_PROPERTY(resname,Resname)
        DEF_PROPERTY(chain,Chain)
        DEF_PROPERTY(name,Name)
        DEF_PROPERTY(mass,Mass)
        DEF_PROPERTY(charge,Charge)
        DEF_PROPERTY(beta,Beta)
        DEF_PROPERTY(occupancy,Occupancy)
        DEF_PROPERTY(resid,Resid)
        DEF_PROPERTY(index,Index)
        DEF_PROPERTY(resindex,Resindex)
        DEF_PROPERTY(tag,Tag)
    ;
}
