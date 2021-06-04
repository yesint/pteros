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





#include "distance_search_contacts_1sel.h"
#include "pteros/core/pteros_error.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;


DistanceSearchContacts1sel::DistanceSearchContacts1sel(float d,
                                                             const Selection& sel,
                                                             std::vector<Eigen::Vector2i>& res_pairs,
                                                             std::vector<float>& res_distances,
                                                             bool absolute_index,
                                                             Vector3i_const_ref pbc)
{
    cutoff = d;
    periodic_dims = pbc;
    is_periodic = (pbc.array()!=0).any();
    abs_index = absolute_index;
    box = sel.box();
    pairs = &res_pairs;
    distances = &res_distances;

    create_grid(sel);

    if(is_periodic){
        grid1.populate_periodic(sel,box,periodic_dims,abs_index);
    } else {
        grid1.populate(sel,min,max,abs_index);
    }

    do_search();
}

void DistanceSearchContacts1sel::search_planned_pair(const PlannedPair& pair,
                                                        std::vector<Vector2i> &pairs_buffer,
                                                        std::vector<float> &distances_buffer)
{
    if(pair.c1==pair.c2){
        // Inside cell
        search_inside_cell(pair,grid1,pairs_buffer,distances_buffer);
    } else {
        // Between cells
        search_between_cells(pair,grid1,grid1,pairs_buffer,distances_buffer);
    }
}




