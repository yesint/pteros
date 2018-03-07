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


#include "distance_search_within_sel.h"
#include "pteros/core/pteros_error.h"
#include "search_utils.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;

Distance_search_within_sel::Distance_search_within_sel(float d,
                            const Selection &src,
                            const Selection &target,
                            std::vector<int>& res,
                            bool include_self,
                            bool periodic){

    cutoff = d;
    is_periodic = periodic;    

    // Get current box
    box = src.box();

    //------------
    // Grid creation part
    //------------

    // Determine bounding box
    if(!is_periodic){
        // Get the minmax of each selection
        Vector3f min1,min2,max1,max2;

        src.minmax(min1,max1);
        target.minmax(min2,max2);
        // Add a "halo: of size cutoff for each of them
        min1.array() -= cutoff;
        max1.array() += cutoff;
        min2.array() -= cutoff;
        max2.array() += cutoff;

        // Find true bounding box
        for(int i=0;i<3;++i){
            overlap_1d(min1(i),max1(i),min2(i),max2(i),min(i),max(i));
            // If no overlap just exit
            if(max(i)==min(i)) return;
        }

    } else {
        // Check if we have periodicity
        if(!box.is_periodic())
            throw Pteros_error("Asked for pbc in within selection, but there is no periodic box!");
        // Set dimensions of the current unit cell
        min.fill(0.0);
        max = box.extents();
    }

    if(src.size()<10 || target.size()<10){
        set_grid_size(min,max, std::max(src.size(),target.size()),box);
    } else {
        set_grid_size(min,max, std::min(src.size(),target.size()),box);
    }

    //set_grid_size(min,max, std::min(src.size(),target.size()),box);

    // Allocate both grids
    grid1.resize(NgridX,NgridY,NgridZ);
    grid2.resize(NgridX,NgridY,NgridZ);

    // Fill grids
    abs_index = false; // Force local indexes!
    if(is_periodic){
        grid1.populate_periodic(src,box,abs_index);
        grid2.populate_periodic(target,box,abs_index);
    } else {
        grid1.populate(src,min,max,abs_index);
        grid2.populate(target,min,max,abs_index);
    }

    do_search(src.size());

    // Here we have to force absolute indexes instead!
    abs_index = true;
    used_to_result(res,include_self,src,target);
}

