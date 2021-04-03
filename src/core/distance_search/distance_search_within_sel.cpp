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


#include "distance_search_within_sel.h"
#include "pteros/core/pteros_error.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;

DistanceSearchWithinSel::DistanceSearchWithinSel(float d,
                            const Selection &src,
                            const Selection &target,
                            std::vector<int>& res,
                            bool include_self,
                            Vector3i_const_ref pbc)
{
    cutoff = d;
    periodic_dims = pbc;
    is_periodic = (pbc.array()!=0).any();
    abs_index = true; // Absolute index is enforced here!
    box = src.box();

    res.clear();

    create_grids(src,target);

    if(is_periodic){
        grid1.populate_periodic(src,box,periodic_dims,abs_index);
        grid2.populate_periodic(target,box,periodic_dims,abs_index);
    } else {
        grid1.populate(src,min,max,abs_index);
        grid2.populate(target,min,max,abs_index);
    }

    do_search();    

    // Elements in set are unique already, need to copy to result and sort
    copy(result.begin(),result.end(),back_inserter(res));
    sort(res.begin(),res.end());

    if(!include_self){
        auto tmp = res;
        res.clear();
        set_difference(tmp.begin(),tmp.end(),target.index_begin(),target.index_end(),back_inserter(res));
    }
}





