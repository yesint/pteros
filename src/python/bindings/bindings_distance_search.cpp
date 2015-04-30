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

#include "pteros/python/bindings_util.h"
#include "bindings_distance_search.h"
#include "pteros/core/distance_search.h"
#include <boost/python/make_constructor.hpp>

using namespace std;
using namespace boost::python;
using namespace pteros;
using namespace Eigen;


boost::shared_ptr<Distance_search_within> distance_search_within_init1(float d,
                                                                      const Selection& src,
                                                                      bool absolute_index,
                                                                      bool periodic)
{
    boost::shared_ptr<Distance_search_within> g(new Distance_search_within(d,src,absolute_index,periodic));
    return g;
}

boost::shared_ptr<Distance_search_within> distance_search_within_init2(float d,
                                                                      const Selection& src,
                                                                      bool absolute_index)
{
    boost::shared_ptr<Distance_search_within> g(new Distance_search_within(d,src,absolute_index,false));
    return g;
}

boost::shared_ptr<Distance_search_within> distance_search_within_init3(float d,
                                                                      const Selection& src)
{
    boost::shared_ptr<Distance_search_within> g(new Distance_search_within(d,src,false,false));
    return g;
}


boost::python::list distance_search_within1(Distance_search_within* g, PyObject* coord){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,c,coord)
    vector<int> r;
    g->search_within(c,r);

    boost::python::list res;
    for(int i=0;i<r.size();++i){
        res.append(r[i]);
    }
    return res;
}

boost::python::list distance_search_within2(Distance_search_within* g,
                                            const Selection& target,
                                            bool include_self){
    vector<int> r;
    g->search_within(target,r,include_self);

    boost::python::list res;
    for(int i=0;i<r.size();++i){
        res.append(r[i]);
    }
    return res;
}

boost::python::list distance_search_within3(Distance_search_within* g,
                                            const Selection& target)
{
    return distance_search_within2(g,target,true);
}

boost::python::tuple search_contacts1(float d,
                     const Selection& sel,
                     bool absolute_index = false,
                     bool periodic = false,
                     bool do_dist = false)
{
    std::vector<Eigen::Vector2i> pairs;
    std::vector<float> dist_vec;
    if(do_dist){
        search_contacts(d,sel,pairs,absolute_index,periodic,&dist_vec);
    } else {
        search_contacts(d,sel,pairs,absolute_index,periodic,nullptr);
    }

    npy_intp dims[2] = {pairs.size(),2};
    PyObject* p1 = PyArray_SimpleNew(2, dims, PyArray_INT);
    Map<MatrixXi>((int*)PyArray_DATA(p1),2,pairs.size()) = Map<MatrixXi>((int*)(pairs.data()),2,pairs.size());

    if(do_dist){
        CREATE_PYARRAY_1D_AND_MAP(p2,VectorXf,v,dist_vec.size())
        v = Map<VectorXf>(dist_vec.data(),dist_vec.size());
        return boost::python::make_tuple(handle<>(p1),handle<>(p2));
    } else {
        return boost::python::make_tuple(handle<>(p1),boost::python::object());
    }
}

BOOST_PYTHON_FUNCTION_OVERLOADS(search_contacts1_overloads, search_contacts1, 2, 5)

boost::python::tuple search_contacts2(float d,
                     const Selection& sel1,
                     const Selection& sel2,
                     bool absolute_index = false,
                     bool periodic = false,
                     bool do_dist = false)
{
    std::vector<Eigen::Vector2i> pairs;
    std::vector<float> dist_vec;
    if(do_dist){
        search_contacts(d,sel1,sel2,pairs,absolute_index,periodic,&dist_vec);
    } else {
        search_contacts(d,sel1,sel2,pairs,absolute_index,periodic,nullptr);
    }

    npy_intp dims[2] = {pairs.size(),2};
    PyObject* p1 = PyArray_SimpleNew(2, dims, PyArray_INT);
    Map<MatrixXi>((int*)PyArray_DATA(p1),2,pairs.size()) = Map<MatrixXi>((int*)(pairs.data()),2,pairs.size());

    if(do_dist){
        CREATE_PYARRAY_1D_AND_MAP(p2,VectorXf,v,dist_vec.size())
        v = Map<VectorXf>(dist_vec.data(),dist_vec.size());
        return boost::python::make_tuple(handle<>(p1),handle<>(p2));
    } else {
        return boost::python::make_tuple(handle<>(p1),boost::python::object());
    }
}

BOOST_PYTHON_FUNCTION_OVERLOADS(search_contacts2_overloads, search_contacts2, 3, 6)

//------------------------------------------------------

void make_bindings_distance_search(){
    import_array();

    class_<Distance_search_within,boost::shared_ptr<Distance_search_within>,boost::noncopyable >("Distance_search_within", init<>())
            .def("__init__",make_constructor(&distance_search_within_init1))
            .def("__init__",make_constructor(&distance_search_within_init2))
            .def("__init__",make_constructor(&distance_search_within_init3))
            .def("setup",&Distance_search_within::setup)
            .def("search_within",&distance_search_within1)
            .def("search_within",&distance_search_within2)
            .def("search_within",&distance_search_within3)
    ;

    def("search_contacts",&search_contacts1,search_contacts1_overloads());
    def("search_contacts",&search_contacts2,search_contacts2_overloads());
}
