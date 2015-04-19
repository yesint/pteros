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
#include "bindings_grid_search.h"
#include "pteros/core/grid_search.h"
#include <boost/python/make_constructor.hpp>
#include <boost/python/slice.hpp>

using namespace std;
using namespace boost::python;
using namespace pteros;
using namespace Eigen;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(assign_to_grid_overloads, assign_to_grid, 2, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_to_custom_grid_overloads, add_to_custom_grid, 1, 3)

/* Here we are forced to work with boost::shared_ptr.
 * We can't use std::shared_ptr since boost.python doesn't wrap it correctly for now...
 */

//---------------- for one selection --------------------

boost::shared_ptr<Grid_searcher> constructor0(float d, Selection& sel,
                                              boost::python::list& bon,
                                              bool absolute_index,
                                              bool periodic,
                                              boost::python::list& dist){

    vector<Vector2i> bond;
    vector<float> dist_vec;
    boost::shared_ptr<Grid_searcher> g(new Grid_searcher(d,sel,bond,absolute_index,periodic,&dist_vec));

    // Convert list of bonds
    del(bon[slice()]);
    for(int i=0;i<bond.size();++i){
        boost::python::list el;
        el.append(bond[i](0));
        el.append(bond[i](1));
        bon.append(el);
    }

    // Convert dist vec
    del(dist[slice()]);
    for(int i=0;i<dist_vec.size();++i){
        dist.append(dist_vec[i]);
    }

    return g;
}

boost::shared_ptr<Grid_searcher> constructor1(float d, Selection& sel,
                                              boost::python::list& bon,
                                              bool absolute_index,
                                              bool periodic
                                              ){

    vector<Vector2i> bond;
    vector<float> dist_vec;
    boost::shared_ptr<Grid_searcher> g(new Grid_searcher(d,sel,bond,absolute_index,periodic,NULL));
    // Convert list of bonds
    del(bon[slice()]);
    for(int i=0;i<bond.size();++i){
        boost::python::list el;
        el.append(bond[i](0));
        el.append(bond[i](1));
        bon.append(el);
    }

    return g;
}

boost::shared_ptr<Grid_searcher> constructor2(float d, Selection& sel,
                                              boost::python::list& bon,
                                              bool absolute_index
                                              ){
    return constructor1(d,sel,bon,absolute_index,false);
}

boost::shared_ptr<Grid_searcher> constructor3(float d, Selection& sel,
                                              boost::python::list& bon
                                              ){
    return constructor1(d,sel,bon,false,false);
}

//---------------- for two selections --------------------

boost::shared_ptr<Grid_searcher> constructor4(float d, Selection& sel1, Selection& sel2,
                                              boost::python::list& bon,
                                              bool absolute_index,
                                              bool periodic,
                                              boost::python::list& dist){

    vector<Vector2i> bond;
    vector<float> dist_vec;
    boost::shared_ptr<Grid_searcher> g(new Grid_searcher(d,sel1,sel2,bond,absolute_index,periodic,&dist_vec));
    // Convert list of bonds
    del(bon[slice()]);
    for(int i=0;i<bond.size();++i){
        boost::python::list el;
        el.append(bond[i](0));
        el.append(bond[i](1));
        bon.append(el);
    }

    // Convert dist vec
    del(dist[slice()]);
    for(int i=0;i<dist_vec.size();++i){
        dist.append(dist_vec[i]);
    }

    return g;
}

boost::shared_ptr<Grid_searcher> constructor5(float d, Selection& sel1, Selection& sel2,
                                              boost::python::list& bon,
                                              bool absolute_index,
                                              bool periodic
                                              ){

    vector<Vector2i> bond;
    vector<float> dist_vec;
    boost::shared_ptr<Grid_searcher> g(new Grid_searcher(d,sel1,sel2,bond,absolute_index,periodic,NULL));
    // Convert list of bonds
    del(bon[slice()]);
    for(int i=0;i<bond.size();++i){
        boost::python::list el;
        el.append(bond[i](0));
        el.append(bond[i](1));
        bon.append(el);
    }

    return g;
}

boost::shared_ptr<Grid_searcher> constructor6(float d, Selection& sel1,Selection& sel2,
                                              boost::python::list& bon,
                                              bool absolute_index
                                              ){
    return constructor5(d,sel1,sel2,bon,absolute_index,false);
}

boost::shared_ptr<Grid_searcher> constructor7(float d, Selection& sel1,Selection& sel2,
                                              boost::python::list& bon
                                              ){
    return constructor5(d,sel1,sel2,bon,false,false);
}

//---------------- for within selections --------------------
void search_within1(Grid_searcher* g, PyObject* coor, boost::python::list& bon){
    MAP_EIGEN_TO_PYTHON_F(Vector3f,c,coor)
    std::vector<int> b;
    Vector3f _c;
    _c = c;
    g->search_within(_c,b);
    del(bon[slice()]);
    for(int i=0;i<b.size();++i){
        bon.append(b[i]);
    }
}

void search_within2(Grid_searcher* g, Selection& target, boost::python::list& bon, bool include_self){
    std::vector<int> b;
    g->search_within(target,b,include_self);
    del(bon[slice()]);
    for(int i=0;i<b.size();++i){
        bon.append(b[i]);
    }
}

void search_within3(Grid_searcher* g, Selection& target, boost::python::list& bon){
    search_within2(g,target,bon,true);
}

/*
boost::python::list cell_of_custom_grid(Grid_searcher* g, int x, int y, int z){
    vector<int> l;
    l = g->cell_of_custom_grid(x,y,z);
    boost::python::list ret;
    for(int i=0;i<l.size();++i){
        ret.append(l[i]);
    }
    return ret;
}
*/

//------------------------------------------------------

void make_bindings_grid_search(){

    class_<Grid_searcher,boost::shared_ptr<Grid_searcher> >("Grid_searcher", init<>())
            .def("__init__",make_constructor(&constructor0))
            .def("__init__",make_constructor(&constructor1))
            .def("__init__",make_constructor(&constructor2))
            .def("__init__",make_constructor(&constructor3))
            .def("__init__",make_constructor(&constructor4))
            .def("__init__",make_constructor(&constructor5))
            .def("__init__",make_constructor(&constructor6))
            .def("__init__",make_constructor(&constructor7))
            .def("assign_to_grid",&Grid_searcher::assign_to_grid,assign_to_grid_overloads())
            .def("search_within",&search_within1)
            .def("search_within",&search_within2)
            .def("search_within",&search_within3)
            .def("create_custom_grid",&Grid_searcher::create_custom_grid)
            .def("clear_custom_grid",&Grid_searcher::clear_custom_grid)
            .def("add_to_custom_grid",&Grid_searcher::add_to_custom_grid,add_to_custom_grid_overloads())
            //.def("cell_of_custom_grid",&cell_of_custom_grid)
    ;
}
