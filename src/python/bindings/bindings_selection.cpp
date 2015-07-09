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

#include "bindings_system.h"
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/python/bindings_util.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace pteros;
using namespace std;
using namespace Eigen;
namespace bp = boost::python;
using namespace boost::python;

/**********************
  Wrappers for Selection
***********************/
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(fit_trajectory_overloads, fit_trajectory, 0, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_average_overloads, get_average, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(principal_orient_overloads, principal_orient, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(non_bond_energy_overloads, non_bond_energy, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(rmsd_overloads, rmsd, 1, 2)


PyObject* Selection_get_xyz(Selection* s){    
    CREATE_PYARRAY_2D_AND_MAP_F(p,MatrixXf,data,3,s->size())
    s->get_xyz(data);
    return boost::python::incref(p);
}

void Selection_set_xyz(Selection* s, PyObject* data){
    MAP_EIGEN_TO_PYTHON_F(MatrixXf,m,data)
    s->set_xyz(m);
}

PyObject* Selection_get_average(Selection* s, int b=0, int e=-1){
    CREATE_PYARRAY_2D_AND_MAP_F(p,MatrixXf,data,3,s->size())
    data = s->average_structure(b,e);
    return boost::python::incref(p);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Selection_get_average_overloads, Selection_get_average, 1, 3)

PyObject* Selection_get_mass(Selection* s){
    CREATE_PYARRAY_1D_AND_MAP_F(p,VectorXf,r,s->size())
    r = Map<VectorXf>(s->get_mass().data(),s->size());
    return incref(p);
}

void Selection_set_mass1(Selection* s, PyObject* data){
    MAP_EIGEN_TO_PYTHON_F(VectorXf,v,data)
    vector<float> r(v.data(),v.data()+v.size());
    s->set_mass(r);
}

void Selection_set_mass2(Selection* s, float data){
    s->set_mass(data);
}

PyObject* Selection_get_beta(Selection* s){
    CREATE_PYARRAY_1D_AND_MAP_F(p,VectorXf,r,s->size())
    r = Map<VectorXf>(s->get_beta().data(),s->size());
    return incref(p);
}

void Selection_set_beta1(Selection* s, PyObject* data){
    MAP_EIGEN_TO_PYTHON_F(VectorXf,v,data)
    vector<float> r(v.data(),v.data()+v.size());
    s->set_beta(r);
}

void Selection_set_beta2(Selection* s, float data){
    s->set_beta(data);
}

PyObject* Selection_get_occupancy(Selection* s){
    CREATE_PYARRAY_1D_AND_MAP_F(p,VectorXf,r,s->size())
    r = Map<VectorXf>(s->get_occupancy().data(),s->size());
    return incref(p);
}

void Selection_set_occupancy1(Selection* s, PyObject* data){
    MAP_EIGEN_TO_PYTHON_F(VectorXf,v,data)
    vector<float> r(v.data(),v.data()+v.size());
    s->set_occupancy(r);
}

void Selection_set_occupancy2(Selection* s, float data){
    s->set_occupancy(data);
}

PyObject* Selection_get_traj(Selection* s, int ind, int b=0, int e=-1){
    CREATE_PYARRAY_2D_AND_MAP_F(p,MatrixXf,data,3,s->get_system()->num_frames())
    data = s->atom_traj(ind,b,e);
    return boost::python::incref(p);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Selection_get_traj_overloads, Selection_get_traj, 2, 4)

boost::python::tuple Selection_minmax(Selection* s){
    CREATE_PYARRAY_1D_AND_MAP_F(pmin,Vector3f,min,3)
    CREATE_PYARRAY_1D_AND_MAP_F(pmax,Vector3f,max,3)
    s->minmax(min,max);
    return boost::python::make_tuple(handle<>(pmin),handle<>(pmax));
}

void Selection_translate(Selection* s, PyObject* vec){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,v,vec)
    s->translate(v);
}

void Selection_rotate_3_arg(Selection* s, PyObject* ar1, PyObject* ar2, PyObject* ar3){
    if(PyArray_Check(ar1)){
        MAP_EIGEN_TO_PYTHON_F(Vector3f,dir,ar1)
        MAP_EIGEN_TO_PYTHON_F(Vector3f,piv,ar3)
        s->rotate(dir,extract<float>(ar2),piv);
    } else {
        MAP_EIGEN_TO_PYTHON_F(Vector3f,piv,ar3)
        s->rotate(extract<int>(ar2),extract<float>(ar2),piv);
    }
}

void Selection_rotate_1_arg(Selection* s, PyObject* ar1){
    MAP_EIGEN_TO_PYTHON_F(Matrix3f,matr,ar1)
    s->rotate(matr);
}

void Selection_rotate_2_arg(Selection* s, PyObject* ar1, PyObject* ar2){
    if(PyArray_Check(ar1)){
        MAP_EIGEN_TO_PYTHON_F(Vector3f,ang,ar1)
        MAP_EIGEN_TO_PYTHON_F(Vector3f,piv,ar2)
        s->rotate(ang,piv);
    } else {
        s->rotate(extract<int>(ar1),extract<float>(ar2));
    }
}


PyObject* fit_transform_py(Selection& sel1, Selection& sel2){
    CREATE_PYARRAY_2D_AND_MAP_F(p,Matrix4f,m,4,4)
    m = fit_transform(sel1,sel2).matrix();
    return boost::python::incref(p);
}

PyObject* Selection_principal_transform(Selection* sel, bool periodic=false){
    CREATE_PYARRAY_2D_AND_MAP_F(p,Matrix4f,m,4,4)
    m = sel->principal_transform(periodic).matrix();
    return boost::python::incref(p);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Selection_principal_transform_overloads, Selection_principal_transform, 1, 2)

void Selection_write_normal(Selection* s, std::string fname, int b=-1,int e=-1){
    s->write(fname,b,e);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(write_overloads, Selection_write_normal, 2, 4)

void Selection_apply_transform(Selection* s, PyObject* t){
    MAP_EIGEN_TO_PYTHON_F(Matrix4f,m,t)
    Affine3f tr(m);
    s->apply_transform(tr);
}

PyObject* Selection_getXYZ1(Selection* s, int ind){
    CREATE_PYARRAY_1D_AND_MAP_F(p,Vector3f,data,3)
    data = s->XYZ(ind);
    return boost::python::incref(p);
}

PyObject* Selection_getXYZ2(Selection* s, int ind, int fr){
    CREATE_PYARRAY_1D_AND_MAP_F(p,Vector3f,data,3)
    data = s->XYZ(ind,fr);
    return boost::python::incref(p);
}

void Selection_setXYZ1(Selection* s, int ind, PyObject* obj){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,data,obj)
    s->XYZ(ind) = data;
}

void Selection_setXYZ2(Selection* s, int ind, int fr, PyObject* obj){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,data,obj)
    s->XYZ(ind,fr) = data;
}

PyObject* Selection_get_index(Selection* s){
    CREATE_PYARRAY_1D_AND_MAP_I(p,VectorXi,r,s->size())
    r = Map<VectorXi>(s->get_index().data(),s->size());
    return incref(p);
}

boost::python::list Selection_get_chain(Selection* s){
    boost::python::list l;    
    for(int i=0;i<s->size();++i) l.append(s->Chain(i));
    return l;
}

void Selection_set_chain1(Selection* s, boost::python::list& data){
    int n = len(data);
    vector<char> r;
    r.resize(n);
    for(int i=0;i<r.size();++i) r[i] = extract<char>(data[i]);
    s->set_chain(r);
}

void Selection_set_chain2(Selection* s, char data){
    s->set_chain(data);
}

boost::python::list Selection_get_unique_chain(Selection* s){
    boost::python::list l;
    vector<char> r = s->get_unique_chain();
    for(int i=0;i<r.size();++i) l.append(r[i]);
    return l;
}

boost::python::list Selection_get_unique_resid(Selection* s){
    boost::python::list l;
    vector<int> r = s->get_unique_resid();
    for(int i=0;i<r.size();++i) l.append(r[i]);
    return l;
}

boost::python::list Selection_get_unique_resindex(Selection* s){
    boost::python::list l;
    vector<int> r = s->get_unique_resindex();
    for(int i=0;i<r.size();++i) l.append(r[i]);
    return l;
}


PyObject* Selection_get_resid(Selection* s){
    CREATE_PYARRAY_1D_AND_MAP_I(p,VectorXi,r,s->size())
    r = Map<VectorXi>(s->get_resid().data(),s->size());
    return incref(p);

}

void Selection_set_resid1(Selection* s, PyObject* data){
    MAP_EIGEN_TO_PYTHON_I(VectorXi,v,data)
    vector<int> r(v.data(),v.data()+v.size());
    s->set_resid(r);
}

void Selection_set_resid2(Selection* s, int data){
    s->set_resid(data);
}

PyObject* Selection_get_resindex(Selection* s){
    CREATE_PYARRAY_1D_AND_MAP_I(p,VectorXi,r,s->size())
    r = Map<VectorXi>(s->get_resindex().data(),s->size());
    return incref(p);
}

boost::python::list Selection_get_name(Selection* s){
    boost::python::list l;    
    for(int i=0;i<s->size();++i) l.append(s->Name(i).c_str());
    return l;
}

void Selection_set_name1(Selection* s, boost::python::list& data){
    int n = len(data);
    vector<string> r;
    r.resize(n);
    for(int i=0;i<r.size();++i) r[i] = extract<string>(data[i]);
    s->set_name(r);
}

void Selection_set_name2(Selection* s, string& data){
    s->set_name(data);
}

boost::python::list Selection_get_resname(Selection* s){
    boost::python::list l;    
    for(int i=0;i<s->size();++i) l.append(s->Resname(i).c_str());
    return l;
}

void Selection_set_resname1(Selection* s, boost::python::list& data){
    int n = len(data);
    vector<string> r;
    r.resize(n);
    for(int i=0;i<r.size();++i) r[i] = extract<string>(data[i]);
    s->set_resname(r);
}

void Selection_set_resname2(Selection* s, string data){
    s->set_resname(data);
}


PyObject* Selection_center(Selection* s, bool mass_weighted, bool periodic){
    CREATE_PYARRAY_1D_AND_MAP_F(p,Vector3f,data,3)
    data = s->center(mass_weighted,periodic);
    return boost::python::incref(p);
}


float rmsd_py(Selection& sel1, int fr1, Selection& sel2, int fr2){
    return rmsd(sel1,fr1,sel2,fr2);
}

void fit_py(Selection& sel1, Selection& sel2){
    fit(sel1,sel2);    
}

Energy_components non_bond_energy_py(const Selection& sel1,
                                     const Selection& sel2,
                                     float cutoff = 0.25,
                                     int fr = -1,
                                     bool periodic = true){
    return non_bond_energy(sel1,sel2,cutoff,fr,periodic);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(non_bond_energy_overloads_free, non_bond_energy_py, 2, 5)

boost::python::tuple Selection_inertia(Selection* s, bool periodic=false){
    CREATE_PYARRAY_1D_AND_MAP_F(m,Vector3f,_m,3)
    CREATE_PYARRAY_2D_AND_MAP_F(a,Matrix3f,_a,3,3)
    s->inertia(_m,_a,periodic);
    return boost::python::make_tuple(handle<>(m),handle<>(a));
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Selection_inertia_overloads, Selection_inertia, 1, 2)

boost::python::list Selection_split_by_connectivity(Selection* s, float d){
    vector<Selection> res;
    s->split_by_connectivity(d,res);
    boost::python::list l;
    for(int i=0;i<res.size();++i) l.append(res[i]);
    return l;
}

boost::python::list Selection_split_by_residue(Selection* s){
    vector<Selection> res;
    s->split_by_residue(res);
    boost::python::list l;
    for(int i=0;i<res.size();++i) l.append(res[i]);
    return l;
}

boost::python::list Selection_split_by_chain(Selection* s){
    vector<Selection> ch;
    s->split_by_chain(ch);
    boost::python::list l;
    for(int i=0;i<ch.size();++i) l.append(ch[i]);
    return l;
}

boost::python::list Selection_split_by_contiguous_index(Selection* s){
    vector<Selection> ch;
    s->split_by_contiguous_index(ch);
    boost::python::list l;
    for(int i=0;i<ch.size();++i) l.append(ch[i]);
    return l;
}

boost::python::list Selection_split_by_contiguous_residue(Selection* s){
    vector<Selection> ch;
    s->split_by_contiguous_residue(ch);
    boost::python::list l;
    for(int i=0;i<ch.size();++i) l.append(ch[i]);
    return l;
}


boost::python::list Selection_each_residue(Selection* s){
    vector<Selection> res;
    s->each_residue(res);
    boost::python::list l;
    for(int i=0;i<res.size();++i) l.append(res[i]);
    return l;
}

float Selection_rmsd1(Selection* s, int fr){
    return s->rmsd(fr);
}

float Selection_rmsd2(Selection* s, int fr1, int fr2){
    return s->rmsd(fr1,fr2);
}

bp::tuple Selection_sasa(Selection* s, float probe_r = 0.14,
                     bool do_total_volume = false,
                     bool do_area_per_atom = false,
                     bool do_volume_per_atom = false)
{
    float total_v,ret;
    std::vector<float> area_per_atom;
    std::vector<float> volume_per_atom;
    bp::list a_list, v_list;
    bp::object tot_obj; //Defaults to None

    float* tot_ptr = do_total_volume ? &total_v : nullptr;
    std::vector<float>* a_ptr = do_area_per_atom ? &area_per_atom : nullptr;
    std::vector<float>* v_ptr = do_volume_per_atom ? &volume_per_atom : nullptr;

    ret = s->sasa(probe_r,tot_ptr,a_ptr,v_ptr);

    if(tot_ptr) tot_obj = bp::object(total_v);
    if(a_ptr) for(int i=0;i<area_per_atom.size();++i) a_list.append(area_per_atom[i]);
    if(v_ptr) for(int i=0;i<volume_per_atom.size();++i) v_list.append(volume_per_atom[i]);

    return boost::python::make_tuple(bp::object(ret),tot_obj,a_list,v_list);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(Selection_sasa_overloads,Selection_sasa,1,5)

boost::shared_ptr<Selection> constructor_from_pyobj0(const System& sys, PyObject* obj){
    cout << "%%% 0" << endl;
    if( PyCallable_Check(obj) ){
        // A lambda, which will call Python function
        auto callback = [&obj](const System& sys, int fr, vector<int>& ind)->void {
            // Call python function
            auto py_ind = bp::call<bp::list>(obj, boost::ref(sys));
            // Convert list to vector
            for(int i=0;i<len(py_ind);++i) ind.push_back( extract<int>(py_ind[i]) );
        };

        boost::shared_ptr<Selection> g(new Selection(sys,callback));
        return g;

    } else if( PyString_Check(obj) ) {
        string str = extract<string>(object( handle<>(borrowed(obj)) ));
        boost::shared_ptr<Selection> g(new Selection(sys,str));
        return g;

    } else if( PySequence_Check(obj) ) {
        bp::object l( handle<>(borrowed(obj)) );
        vector<int> v(len(l));
        for(int i=0;i<len(l);++i) v[i] = extract<int>(l[i]);
        boost::shared_ptr<Selection> g(new Selection(sys,v));
        return g;
    } else {
        throw Pteros_error("Invalid arguments for Selection constructor!");
    }
}

boost::shared_ptr<Selection> constructor_from_pyobj1(const System& sys, PyObject* obj, int fr){
    cout << "%%% 1" << endl;
    if( PyCallable_Check(obj) ){
        // A lambda, which will call Python function
        auto callback = [&obj](const System& sys, int fr, vector<int>& ind)->void {
            // Call python function
            auto py_ind = bp::call<bp::list>(obj, boost::ref(sys), fr);
            // Convert list to vector
            for(int i=0;i<len(py_ind);++i) ind.push_back( extract<int>(py_ind[i]) );
        };

        boost::shared_ptr<Selection> g(new Selection(sys,callback,fr));
        return g;

    } else if( PyString_Check(obj) ) {
        string str = extract<string>(object( handle<>(borrowed(obj)) ));        
        boost::shared_ptr<Selection> g(new Selection(sys,str,fr));
        return g;

    } else if( PySequence_Check(obj) ) {
        throw Pteros_error("Specified target frame for selection from a list of indexes! This doesn't make sense.");
        bp::object l( handle<>(borrowed(obj)) );
        vector<int> v(len(l));
        for(int i=0;i<len(l);++i) v[i] = extract<int>(l[i]);
        boost::shared_ptr<Selection> g(new Selection(sys,v));
        return g;

    } else if( PyInt_Check(obj) ) {
        int i = extract<int>(object( handle<>(borrowed(obj)) ));
        boost::shared_ptr<Selection> g(new Selection(sys,i,fr));
        return g;

    } else {
        throw Pteros_error("Invalid arguments for Selection constructor!");
    }
}


void Selection_modify_list1(Selection* sel, const bp::list& l){
    vector<int> v(len(l));
    for(int i=0;i<len(l);++i) v[i] = extract<int>(l[i]);
    sel->modify(v);
}

void Selection_modify_list2(Selection* sel, const System& sys, const bp::list& l){
    vector<int> v(len(l));
    for(int i=0;i<len(l);++i) v[i] = extract<int>(l[i]);
    sel->modify(sys,v);
}

void Selection_modify_from_pyobj3(Selection* sel, const System& sys, PyObject* obj, int fr){
    if( PyCallable_Check(obj) ){
        // A lambda, which will call Python function
        auto callback = [&obj](const System& sys, int fr, vector<int>& ind)->void {
            // Call python function
            auto py_ind = bp::call<bp::list>(obj, boost::ref(sys), fr);
            // Convert list to vector
            for(int i=0;i<len(py_ind);++i) ind.push_back( extract<int>(py_ind[i]) );
        };
        sel->modify(sys,callback,fr);

    } else if( PyString_Check(obj) ) {
        string str = extract<string>(object( handle<>(borrowed(obj)) ));
        sel->modify(sys,str,fr);

    } else if( PySequence_Check(obj) ) {
        throw Pteros_error("Specified target frame for selection from a list of indexes! This doesn't make sense.");
        bp::object l( handle<>(borrowed(obj)) );
        vector<int> v(len(l));
        for(int i=0;i<len(l);++i) v[i] = extract<int>(l[i]);
        sel->modify(sys,v);
    } else if( PyInt_Check(obj) ) {
        int i = extract<int>(object( handle<>(borrowed(obj)) ));
        sel->modify(sys,i,fr);
    } else {
        throw Pteros_error("Wrong arguments for selection!");
    }
}


void Selection_modify_from_pyobj2(Selection* sel, const System& sys, PyObject* obj){        
    if( PyCallable_Check(obj) ){
        // A lambda, which will call Python function
        auto callback = [&obj](const System& sys, int fr, vector<int>& ind)->void {
            // Call python function
            auto py_ind = bp::call<bp::list>(obj, boost::ref(sys), fr);
            // Convert list to vector
            for(int i=0;i<len(py_ind);++i) ind.push_back( extract<int>(py_ind[i]) );
        };
        sel->modify(sys,callback);

    } else if( PyString_Check(obj) ) {
        string str = extract<string>(object( handle<>(borrowed(obj)) ));        
        sel->modify(sys,str);

    } else if( PySequence_Check(obj) ) {
        bp::object l( handle<>(borrowed(obj)) );
        vector<int> v(len(l));
        for(int i=0;i<len(l);++i) v[i] = extract<int>(l[i]);
        sel->modify(sys,v);
    } else {
        throw Pteros_error("Wrong arguments for selection!");
    }
}

void Selection_modify_from_pyobj2_def(Selection* sel, PyObject* obj){
    Selection_modify_from_pyobj2(sel,*sel->get_system(),obj);
}

void Selection_modify_from_pyobj3_def(Selection* sel, PyObject* obj, int fr){
    Selection_modify_from_pyobj3(sel,*sel->get_system(),obj,fr);
}


float Selection_distance1(Selection* s, int i, int j, bool pbc, PyObject* dims){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims)
    return s->distance(i,j,pbc,dim);
}

float Selection_distance2(Selection* s, int i, int j, bool pbc){
    return s->distance(i,j,pbc);
}

float Selection_distance3(Selection* s, int i, int j, int fr){
    return s->distance(i,j);
}

float Selection_angle1(Selection* s, int i, int j, int k, bool pbc, PyObject* dims){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims)
    return s->angle(i,j,k,pbc,dim);
}

float Selection_angle2(Selection* s, int i, int j, int k, bool pbc){
    return s->angle(i,j,k,pbc);
}

float Selection_angle3(Selection* s, int i, int j, int k, int fr){
    return s->angle(i,j,k,fr);
}

float Selection_dihedral1(Selection* s, int i, int j, int k, int l, bool pbc, PyObject* dims){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims)
    return s->dihedral(i,j,k,l,pbc,dim);
}

float Selection_dihedral2(Selection* s, int i, int j, int k, int l, bool pbc){
    return s->dihedral(i,j,k,l,pbc);
}

float Selection_dihedral3(Selection* s, int i, int j, int k, int l){
    return s->dihedral(i,j,k,l);
}

void Selection_wrap0(Selection* s){
    s->wrap();
}

void Selection_wrap1(Selection* s, PyObject* dims){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims)
    s->wrap(dim);
}

void Selection_unwrap0(Selection* s){
    s->unwrap();
}

void Selection_unwrap1(Selection* s, PyObject* dims){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims)
    s->unwrap(dim);
}

void Selection_unwrap_bonds0(Selection* s){
    s->unwrap_bonds();
}

void Selection_unwrap_bonds1(Selection* s, float d){
    s->unwrap_bonds(d);
}

void Selection_unwrap_bonds2(Selection* s, float d, int leading_index){
    s->unwrap_bonds(d,leading_index);
}

void Selection_unwrap_bonds3(Selection* s, float d, int leading_index, PyObject* dims){
    MAP_EIGEN_TO_PYTHON_I(Vector3i,dim,dims)
    s->unwrap_bonds(d,leading_index,dim);
}

//-------------------------------------------------------

// Macros to wrap an inline accessor function
#define WRAP_ACCESSOR_FR(_out_type, _func) \
    _out_type Selection_get##_func##1(Selection* s, int arg){ return s->_func(arg); } \
    void Selection_set##_func##1(Selection* s, int arg, _out_type val){ s->_func(arg) = val; } \
    _out_type Selection_get##_func##2(Selection* s, int arg, int fr){ return s->_func(arg,fr); } \
    void Selection_set##_func##2(Selection* s, int arg, int fr, _out_type val){ s->_func(arg,fr) = val; }

#define WRAP_ACCESSOR(_out_type, _func) \
    _out_type Selection_get##_func(Selection* s, int arg){ return s->_func(arg); } \
    void Selection_set##_func(Selection* s, int arg, _out_type val){ s->_func(arg) = val; }

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
WRAP_ACCESSOR(int,Index)
WRAP_ACCESSOR(std::string,Tag)
WRAP_ACCESSOR(int,Resindex)

Periodic_box Selection_getBox(Selection* sel){
    return sel->Box();
}

void Selection_setBox(Selection* sel, const Periodic_box& b){
    sel->Box() = b;
}

float Selection_getTime(Selection* sel){
    return sel->Time();
}

void Selection_setTime(Selection* sel, const float& t){
    sel->Time() = t;
}

// Iteration support
// For some unknown reason direct wrapping of Selection::iterator does not work properly
// So we create custom iterator class, which satisfyes Python iterator protocol
// and expose it
struct Selection_iter {
    Selection_iter(Selection& sel): it(sel.begin()), it_end(sel.end()) {}

    Atom_proxy next(){
        Selection::iterator i = it;
        if(i==it_end){
            PyErr_SetString(PyExc_StopIteration, "No more data");
            boost::python::throw_error_already_set();
        }
        it++;
        return *i;
    }

    Selection::iterator it, it_end;
};

// Returns an iterator object from __iter__
Selection_iter Selection_get_iter(Selection& sel){
    return Selection_iter(sel);
}

//-------------------
// Indexing

Atom_proxy Selection_getitem(Selection* sel, int i){
    return (*sel)[i];
}



void make_bindings_Selection(){
    import_array();

    class_<Selection_iter>("_Selection_iter", no_init)
        .def("next",&Selection_iter::next)
    ;

    def("rmsd",&rmsd_py);
    def("fit",&fit_py);
    def("fit_transform",&fit_transform_py);
    def("non_bond_energy",&non_bond_energy_py,non_bond_energy_overloads_free());

    class_<Selection>("Selection", init<>())        
        .def(init<System&>() )
        .def(init<const Selection&>() )

        .def("__init__",make_constructor(&constructor_from_pyobj0))
        .def("__init__",make_constructor(&constructor_from_pyobj1))

        .def("size",&Selection::size)

        // Modification of existing selection

        .def("append", static_cast<void(Selection::*)   (const Selection&)      >(&Selection::append) )
        .def("append", static_cast<void(Selection::*)   (int)                   >(&Selection::append) )        

        .def("remove", static_cast<void(Selection::*)   (const Selection&)      >(&Selection::remove) )
        .def("remove", static_cast<void(Selection::*)   (int)                   >(&Selection::remove) )

        .def("set_system", &Selection::set_system)

        .def("modify", &Selection_modify_from_pyobj2)
        .def("modify", &Selection_modify_from_pyobj2_def)
        .def("modify", &Selection_modify_from_pyobj3)
        .def("modify", &Selection_modify_from_pyobj3_def)

        .def("apply",&Selection::apply)
        .def("update",&Selection::update)

        .def("get_frame",&Selection::get_frame)
        .def("set_frame",&Selection::set_frame)

        .def("clear",&Selection::clear)        

        .def("get_system",&Selection::get_system, return_value_policy<reference_existing_object>())

        .def("get_text",&Selection::get_text)

        .def("get_index",&Selection_get_index)

        .def("get_chain",&Selection_get_chain)
        .def("set_chain",&Selection_set_chain1)
        .def("set_chain",&Selection_set_chain2)

        .def("get_unique_chain",&Selection_get_unique_chain)

        .def("get_resid",&Selection_get_resid)
        .def("set_resid",&Selection_set_resid1)
        .def("set_resid",&Selection_set_resid2)

        .def("get_unique_resid",&Selection_get_unique_resid)        

        .def("get_resindex",&Selection_get_resindex)                    
        .def("get_unique_resindex",&Selection_get_unique_resindex)

        .def("get_name",&Selection_get_name)
        .def("set_name",&Selection_set_name1)
        .def("set_name",&Selection_set_name2)

        .def("get_resname",&Selection_get_resname)
        .def("set_resname",&Selection_set_resname1)
        .def("set_resname",&Selection_set_resname2)

        .def("get_xyz",&Selection_get_xyz)
        .def("set_xyz",&Selection_set_xyz)

        .def("get_average",&Selection_get_average, Selection_get_average_overloads())

        .def("get_mass",&Selection_get_mass)
        .def("set_mass",&Selection_set_mass1)
        .def("set_mass",&Selection_set_mass2)

        .def("get_beta",&Selection_get_beta)
        .def("set_beta",&Selection_set_beta1)
        .def("set_beta",&Selection_set_beta2)

        .def("get_occupancy",&Selection_get_occupancy)
        .def("set_occupancy",&Selection_set_occupancy1)
        .def("set_occupancy",&Selection_set_occupancy2)

        .def("get_traj",&Selection_get_traj, Selection_get_traj_overloads())

        .def("center",&Selection_center, (bp::arg("mass_weighted")=false,bp::arg("periodic")=false) )

        .def("minmax",&Selection_minmax)

        .def("distance",&Selection_distance1)
        .def("distance",&Selection_distance2)
        .def("distance",&Selection_distance3)

        .def("angle",&Selection_angle1)
        .def("angle",&Selection_angle2)
        .def("angle",&Selection_angle3)

        .def("dihedral",&Selection_dihedral1)
        .def("dihedral",&Selection_dihedral2)
        .def("dihedral",&Selection_dihedral3)

        .def("translate",&Selection_translate)

        .def("rotate",&Selection_rotate_1_arg)
        .def("rotate",&Selection_rotate_2_arg)
        .def("rotate",&Selection_rotate_3_arg)

        .def("rmsd", &Selection_rmsd1 )
        .def("rmsd", &Selection_rmsd2 )

        .def("fit_trajectory",&Selection::fit_trajectory, fit_trajectory_overloads())
        .def("apply_transform",&Selection_apply_transform)

        .def("write",&Selection_write_normal, write_overloads())

        // inertia return a tuple of (moments,axes)
        .def("inertia",&Selection_inertia, Selection_inertia_overloads())

        .def("principal_transform",&Selection_principal_transform, Selection_principal_transform_overloads())

        .def("principal_orient",&Selection::principal_orient, principal_orient_overloads())

        .def("wrap",&Selection_wrap0)
        .def("wrap",&Selection_wrap1)
        .def("unwrap",&Selection_unwrap0)
        .def("unwrap",&Selection_unwrap1)
        .def("unwrap_bonds",&Selection_unwrap_bonds0)
        .def("unwrap_bonds",&Selection_unwrap_bonds1)
        .def("unwrap_bonds",&Selection_unwrap_bonds2)
        .def("unwrap_bonds",&Selection_unwrap_bonds3)

        .def("text_based",&Selection::text_based)
        .def("coord_dependent",&Selection::coord_dependent)
        .def("flatten",&Selection::flatten)

        .def("non_bond_energy",&Selection::non_bond_energy, non_bond_energy_overloads())

        // Splitting functions
        .def("split_by_connectivity",&Selection_split_by_connectivity)
        .def("split_by_residue",&Selection_split_by_residue)
        .def("split_by_chain",&Selection_split_by_chain)
        .def("split_by_contiguous_index",&Selection_split_by_contiguous_index)
        .def("split_by_contiguous_residue",&Selection_split_by_contiguous_residue)
        .def("each_residue",&Selection_each_residue)

        .def("sasa",&Selection_sasa,Selection_sasa_overloads())

        // For coordinate accessors we should use setX instead of just X in Python
        // This is because Python don't respect void in return - all functions
        // with equal number of argumets are considered equivalent thus
        // "float X(int,int)" and "void X(int,float)" become the same function...
        // For non-coordinate accessors this is not needed but used to be consistent
        .def("getX",&Selection_getX1)
        .def("getX",&Selection_getX2)
        .def("setX",&Selection_setX1)
        .def("setX",&Selection_setX2)

        .def("getY",&Selection_getY1)
        .def("getY",&Selection_getY2)
        .def("setY",&Selection_setY1)
        .def("setY",&Selection_setY2)

        .def("getZ",&Selection_getZ1)
        .def("getZ",&Selection_getZ2)
        .def("setZ",&Selection_setZ1)
        .def("setZ",&Selection_setZ2)

        .def("getXYZ",&Selection_getXYZ1)
        .def("getXYZ",&Selection_getXYZ2)
        .def("setXYZ",&Selection_setXYZ1)
        .def("setXYZ",&Selection_setXYZ2)

        //.def("getXYZ1",static_cast<const Eigen::Vector3f&(Selection::*)(int)const>(&Selection::XYZ),return_value_policy<reference_existing_object>())

        .def("getType",&Selection_getType)
        .def("setType",&Selection_setType)

        .def("getType_name",&Selection_getType_name)
        .def("setType_name",&Selection_setType_name)

        .def("getResname",&Selection_getResname)
        .def("setResname",&Selection_setResname)

        .def("getChain",&Selection_getChain)
        .def("setChain",&Selection_setChain)

        .def("getName",&Selection_getName)
        .def("setName",&Selection_setName)

        .def("getMass",&Selection_getMass)
        .def("setMass",&Selection_setMass)

        .def("getCharge",&Selection_getCharge)
        .def("setCharge",&Selection_setCharge)

        .def("getBeta",&Selection_getBeta)
        .def("setBeta",&Selection_setBeta)

        .def("getOccupancy",&Selection_getOccupancy)
        .def("setOccupancy",&Selection_setOccupancy)

        .def("getResid",&Selection_getResid)
        .def("setResid",&Selection_setResid)

        .def("getIndex",&Selection_getIndex)
        .def("setIndex",&Selection_setIndex)

        .def("getResindex",&Selection_getResindex)
        .def("setResindex",&Selection_setResindex)

        .def("getTag",&Selection_getTag)
        .def("setTag",&Selection_setTag)

        .def("getBox",&Selection_getBox)
        .def("setBox",&Selection_setBox)        

        .def("getTime",&Selection_getTime)
        .def("setTime",&Selection_setTime)

        // Iteration protocol support
        .def("__iter__", &Selection_get_iter)
        // Indexing support
        .def("__len__",&Selection::size)
        .def("__getitem__",&Selection_getitem)

        // String conversion
        .def(self_ns::str(self_ns::self))

        // Operators
        .def(self == self)
        .def(self != self)
        .def(self | self)
        .def(self & self)
        .def(self - self)
        .def(~self)

    ;
}
