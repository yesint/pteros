/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#include "pteros/core/distance_search.h"
#include "bindings_util.h"

namespace py = pybind11;
using namespace pteros;
using namespace std;
using namespace Eigen;
using namespace pybind11::literals;


void make_bindings_Distance_search(py::module& m){

    m.def("search_contacts",[](float d,
          const Selection& sel,
          bool absolute_index = false,
          Array3i_const_ref pbc = noPBC )
    {
        std::vector<float>* dist_vec_ptr = new std::vector<float>;
        std::vector<Vector2i>* pairs_ptr = new std::vector<Vector2i>;

        search_contacts(d,sel,*pairs_ptr,*dist_vec_ptr,absolute_index,pbc);

        if(pairs_ptr->size()){
            // Interpret pairs array as 1D array of ints first and convert to py::array
            // Pass size*2 explicitly to ensure correct size
            py::array m = vector_to_array<int>(reinterpret_cast<std::vector<int>*>(pairs_ptr),2*pairs_ptr->size());
            // Reshape into 2D array (no reallocation)
            m.resize(vector<size_t>{pairs_ptr->size(),2});

            // Make py::array for distances
            py::array v = vector_to_array<float>(dist_vec_ptr);
            return py::make_tuple(m,v);
        } else {
            // Empty
            return py::make_tuple(py::list(),py::list());
        }

    }, "d"_a, "sel"_a, "abs_ind"_a=false, "pbc"_a=noPBC);


    m.def("search_contacts",[](float d,
          const Selection& sel1,
          const Selection& sel2,
          bool absolute_index = false,
          Vector3i_const_ref pbc = noPBC )
    {
        std::vector<float>* dist_vec_ptr = new std::vector<float>;
        std::vector<Vector2i>* pairs_ptr = new std::vector<Vector2i>;
        search_contacts(d,sel1,sel2,*pairs_ptr,*dist_vec_ptr,absolute_index,pbc);

        if(pairs_ptr->size()){
            // Interpret pairs array as 1D array of ints first
            // Pass size*2 explicitly to ensure correct size
            py::array m = vector_to_array<int>(reinterpret_cast<std::vector<int>*>(pairs_ptr),2*pairs_ptr->size());
            // Reshape into 2D array (no reallocation)
            m.resize(vector<size_t>{pairs_ptr->size(),2});

            // Make py::array for distances
            py::array v = vector_to_array<float>(dist_vec_ptr);
            return py::make_tuple(m,v);

        } else {
            // Empty
            return py::make_tuple(py::list(),py::list());
        }

    }, "d"_a, "sel1"_a, "sel2"_a, "abs_ind"_a=false, "pbc"_a=noPBC);


    m.def("search_within",[](float d,
          const Selection& src,
          const Selection& target,
          bool include_self = true,
          Vector3i_const_ref pbc = noPBC )
    {
        std::vector<int>* res_ptr = new std::vector<int>;
        search_within(d,src,target,*res_ptr,include_self,pbc);
        return vector_to_array<int>(res_ptr);
    }, "d"_a, "src"_a, "target"_a, "include_self"_a=true, "pbc"_a=noPBC);


    py::class_<DistanceSearchWithin>(m, "DistanceSearchWithin")
            .def(py::init<float,const Selection&,bool,Vector3i_const_ref>(), "d"_a,"src"_a,"abs_ind"_a=false,"pbc"_a=noPBC)
            .def("setup",&DistanceSearchWithin::setup, "d"_a,"src"_a,"abs_ind"_a=false,"pbc"_a=noPBC)

            .def("search_within",[](DistanceSearchWithin* obj, Vector3f_const_ref coord)
                {
                   std::vector<int>* res_ptr = new std::vector<int>;
                   obj->search_within(coord,*res_ptr);
                   return vector_to_array<int>(res_ptr);
                })

            .def("search_within",[](DistanceSearchWithin* obj, const Selection& target, bool include_self)
                {
                    std::vector<int>* res_ptr = new std::vector<int>;
                    obj->search_within(target,*res_ptr,include_self);
                    return vector_to_array<int>(res_ptr);
                },"target"_a, "include_self"_a=true)
    ;
}




