/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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


#include "pteros/core/distance_search_within.h"
#include "pteros/core/pteros_error.h"
#include "distance_search_within_base.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;


class Distance_search_within::Distance_search_within_impl: public Distance_search_within_base {
public:
    Distance_search_within_impl(){}

    /// Constructor for two-stage within searching.
    /// Sets source selection.
    Distance_search_within_impl(float d,
                           const Selection& src,
                           bool absolute_index = false,
                           bool periodic = false)
    {
        setup(d,src,absolute_index,periodic);
    }

    ~Distance_search_within_impl(){}


    void setup(float d,
               const Selection& src,
               bool absolute_index = false,
               bool periodic = false)
    {
        cutoff = d;
        is_periodic = periodic;
        abs_index = absolute_index;
        box = src.box();

        // Grid creation
        if(!is_periodic){
            // Get the minmax of selection
            src.minmax(min,max);
            // Add a "halo: of size cutoff
            min.array() -= cutoff;
            max.array() += cutoff;
        } else {
            // Check if we have periodicity
            if(!box.is_periodic())
                throw Pteros_error("Asked for pbc in within selection, but there is no periodic box!");
            // Set dimensions of the current unit cell
            min.fill(0.0);
            max = box.extents();
        }

        set_grid_size(min,max, src.size(), box);
        // Allocate first grid
        grid1.resize(NgridX,NgridY,NgridZ);

        // Populate first grid
        // We have to force local indexes here in order to allow atomic array to work correctly
        bool saved_abs = abs_index;
        abs_index = false;
        if(is_periodic){
            grid1.populate_periodic(src,box,abs_index);
        } else {
            grid1.populate(src,min,max,abs_index);
        }
        // Restore abs_index
        abs_index = saved_abs;

        src_ptr = const_cast<Selection*>(&src);
    }


    /// Search atoms from source within given distance from given point in space
    void search_within(Vector3f_const_ref coord,
                       std::vector<int> &res)
    {

        // grid1 has parent selection, which is "src". Our point is "target"
        System tmp; // tmp system with just one atom with coordinates coord
        vector<Vector3f> crd{coord};
        vector<Atom> atm(1);
        tmp.atoms_add(atm,crd);
        auto target = tmp.select_all();
        // Allocate second grid of the same size
        grid2.resize(NgridX,NgridY,NgridZ);

        // We have to force local indexes here in order to allow atomic array to work correctly
        bool saved_abs = abs_index;
        abs_index = false;
        if(is_periodic){
            grid2.populate_periodic(target,box,abs_index);
        } else {
            grid2.populate(target,min,max,abs_index);
        }
        // Restore
        abs_index = saved_abs;

        // Now search
        do_search(src_ptr->size());

        //second src_ptr is not used inside. Passed just to satisfy signature
        used_to_result(res,true,*src_ptr,*src_ptr);
    }


    /// Search atoms from source within given distance from target selection
    /// \warning Target must be the subset of source to give meaningful results!
    void search_within(const Selection& target,
                       std::vector<int> &res,
                       bool include_self=true)
    {
        // Allocate second grid of the same size
        grid2.resize(NgridX,NgridY,NgridZ);

        if(is_periodic){
            grid2.populate_periodic(target,box,abs_index);
        } else {
            grid2.populate(target,min,max,abs_index);
        }

        do_search(src_ptr->size());

        used_to_result(res,include_self,*src_ptr,target);
    }

};



Distance_search_within::Distance_search_within()
{
    p = unique_ptr<Distance_search_within_impl>(new Distance_search_within_impl());
}

Distance_search_within::Distance_search_within(float d, const Selection &src, bool absolute_index, bool periodic)
{
    p = unique_ptr<Distance_search_within_impl>(new Distance_search_within_impl(d,src,absolute_index,periodic));
}

Distance_search_within::~Distance_search_within()
{

}

void Distance_search_within::setup(float d, const Selection &src, bool absolute_index, bool periodic)
{
    p->setup(d,src,absolute_index,periodic);
}

void Distance_search_within::search_within(Vector3f_const_ref coord, std::vector<int> &res)
{
    p->search_within(coord,res);
}

void Distance_search_within::search_within(const Selection &target, std::vector<int> &res, bool include_self)
{
    p->search_within(target,res,include_self);
}

