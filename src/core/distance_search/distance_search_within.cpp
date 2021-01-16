/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2020, Semen Yesylevskyy
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
#include "distance_search_within_base.h"

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
                           Vector3i_const_ref pbc = fullPBC)
    {
        setup(d,src,absolute_index,pbc);
    }

    ~Distance_search_within_impl(){}


    void setup(float d,
               const Selection& src,
               bool absolute_index = false,
               Vector3i_const_ref pbc = fullPBC)
    {
        cutoff = d;
        periodic_dims = pbc;
        is_periodic = (pbc.array()!=0).any();
        abs_index = absolute_index;
        box = src.box();

        create_grid(src);

        if(is_periodic){
            grid1.populate_periodic(src,box,periodic_dims,abs_index);
        } else {
            grid1.populate(src,min,max,abs_index);
        }
    }


    /// Search atoms from source within given distance from given point in space
    void search_within(Vector3f_const_ref coord,
                       std::vector<int> &res)
    {
        res.clear();
        // grid1 has parent selection, which is "src". Our point is "target"
        System tmp; // tmp system with just one atom with coordinates coord
        vector<Vector3f> crd{coord};
        vector<Atom> atm(1);
        tmp.atoms_add(atm,crd);
        auto target = tmp.select_all(); // This is our target selection

        // Allocate second grid of the same size
        grid2.resize(Ngrid(0),Ngrid(1),Ngrid(2));

        // We have to force local indexes here in order to allow atomic array to work correctly
        if(is_periodic){
            grid2.populate_periodic(target,box,periodic_dims,abs_index);
        } else {
            grid2.populate(target,min,max,abs_index);
        }

        // Now search
        do_search();

        copy(result.begin(),result.end(),back_inserter(res));
        sort(res.begin(),res.end());
    }


    /// Find atoms from source within given distance from target selection
    void search_within(const Selection& target,
                       std::vector<int> &res,
                       bool include_self=true)
    {
        res.clear();
        // Allocate second grid of the same size
        grid2.resize(Ngrid(0),Ngrid(1),Ngrid(2));

        if(is_periodic){
            grid2.populate_periodic(target,box,periodic_dims,abs_index);
        } else {
            grid2.populate(target,min,max,abs_index);
        }

        do_search();

        if(include_self) result.insert(target.index_begin(),target.index_end());
        // Elements in set are unique already, need to copy to result and sort
        copy(result.begin(),result.end(),back_inserter(res));
        sort(res.begin(),res.end());
    }

};



Distance_search_within::Distance_search_within()
{
    p = unique_ptr<Distance_search_within_impl>(new Distance_search_within_impl());
}

Distance_search_within::Distance_search_within(float d, const Selection &src, bool absolute_index, Vector3i_const_ref pbc)
{
    p = unique_ptr<Distance_search_within_impl>(new Distance_search_within_impl(d,src,absolute_index,pbc));
}

Distance_search_within::~Distance_search_within()
{

}

void Distance_search_within::setup(float d, const Selection &src, bool absolute_index, Vector3i_const_ref pbc)
{
    p->setup(d,src,absolute_index,pbc);
}

void Distance_search_within::search_within(Vector3f_const_ref coord, std::vector<int> &res)
{
    p->search_within(coord,res);
}

void Distance_search_within::search_within(const Selection &target, std::vector<int> &res, bool include_self)
{
    p->search_within(target,res,include_self);
}


