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





#include "distance_search_contacts_2sel.h"
#include "pteros/core/pteros_error.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;


DistanceSearchContacts2sel::DistanceSearchContacts2sel(float d,
                                                             const Selection &sel1,
                                                             const Selection &sel2,
                                                             std::vector<Eigen::Vector2i>& res_pairs,
                                                             std::vector<float> &res_distances,
                                                             bool absolute_index,
                                                             Vector3i_const_ref pbc)
{
    cutoff = d;
    periodic_dims = pbc;
    is_periodic = (pbc.array()!=0).any();
    abs_index = absolute_index;

    if(sel1.get_system() != sel2.get_system())
        throw PterosError("Selections for distance search should be from the same system!");

    box = sel1.box();
    pairs = &res_pairs;
    distances = &res_distances;

    create_grids(sel1,sel2);

    if(is_periodic){
        grid1.populate_periodic(sel1,box,periodic_dims,abs_index);
        grid2.populate_periodic(sel2,box,periodic_dims,abs_index);
    } else {
        grid1.populate(sel1,min,max,abs_index);
        grid2.populate(sel2,min,max,abs_index);
    }

    do_search();
}

void DistanceSearchContacts2sel::search_planned_pair(const PlannedPair &pair,
                                                        std::vector<Vector2i> &pairs_buffer,
                                                        std::vector<float> &distances_buffer)
{
    search_between_cells(pair,grid1,grid2,pairs_buffer,distances_buffer);
    if(pair.c1!=pair.c2) search_between_cells(pair,grid2,grid1,pairs_buffer,distances_buffer);
}





