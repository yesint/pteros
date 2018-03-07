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


#include "distance_search_contacts_2sel.h"
#include "search_utils.h"
#include "pteros/core/pteros_error.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;


Distance_search_contacts_2sel::Distance_search_contacts_2sel(float d,
                                                             const Selection &sel1,
                                                             const Selection &sel2,
                                                             std::vector<Eigen::Vector2i>& pairs,
                                                             bool absolute_index,
                                                             bool periodic,
                                                             std::vector<float> *dist_vec)
{
    result_pairs = &pairs;
    result_distances = dist_vec;
    cutoff = d;
    is_periodic = periodic;
    abs_index = absolute_index;
    box = sel1.box();

    create_grids(sel1,sel2);

    if(is_periodic){
        grid1.populate_periodic(sel1,box,abs_index);
        grid2.populate_periodic(sel2,box,abs_index);
    } else {
        grid1.populate(sel1,min,max,abs_index);
        grid2.populate(sel2,min,max,abs_index);
    }

    do_search();
}

void Distance_search_contacts_2sel::create_grids(const Selection &sel1, const Selection &sel2)
{
    if(!is_periodic){
        // Get the minmax of each selection
        Vector3f min1,min2,max1,max2;

        sel1.minmax(min1,max1);
        sel2.minmax(min2,max2);
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

    set_grid_size(min,max, sel1.size()+sel2.size(), box);
    // Allocate both grids
    grid1.resize(NgridX,NgridY,NgridZ);
    grid2.resize(NgridX,NgridY,NgridZ);
    // Allocate visited array
    visited.resize( boost::extents[NgridX][NgridY][NgridZ] );
}

void Distance_search_contacts_2sel::do_part(int dim, int _b, int _e, std::deque<Vector2i> &bon, std::deque<float> *dist_vec)
{
    Vector3i b(0,0,0);
    Vector3i e(NgridX,NgridY,NgridZ);
    int dim_max = e(dim);
    b(dim)= _b;
    e(dim)= _e;
    int i,j,k,i1,nlist_size;
    Nlist_t nlist; // Local nlist

    int s1,s2,s3;

    for(i=b(0);i<e(0);++i){
        for(j=b(1);j<e(1);++j){
            for(k=b(2);k<e(2);++k){
                // Search in central cell
                // Central cell is always non-periodic
                search_in_pair_of_cells(i,j,k, i,j,k,
                                        grid1,grid2,
                                        bon,dist_vec,
                                        false);
                visited[i][j][k] = true;
                // Get neighbour list locally
                get_nlist(i,j,k,nlist);
                nlist_size = nlist.data.size();
                // Search between this and neighbouring cells
                for(i1=0;i1<nlist_size;++i1){
                    const Vector3i& cell = nlist.data[i1];

                    // If the neighbour is "at left" from the boundary of this part,
                    // ignore it. Only consider dim dimension.
                    if(cell(dim)<b(dim)){
                        continue;
                    }

                    s1 = cell(0);
                    s2 = cell(1);
                    s3 = cell(2);

                    // We only check for visited cells inside local part, not in the "halo"
                    if(    cell(dim)>=b(dim)
                        && cell(dim)<e(dim) ){
                        // cell is inside the local partition
                        if( !visited[s1][s2][s3] ){
                            search_in_pair_of_cells(i,j,k,
                                                    s1,s2,s3,
                                                    grid1, grid2,
                                                    bon,dist_vec,
                                                    nlist.wrapped[i1] && is_periodic);
                            search_in_pair_of_cells(s1,s2,s3,
                                                    i,j,k,
                                                    grid1, grid2,
                                                    bon,dist_vec,
                                                    nlist.wrapped[i1] && is_periodic);
                        }
                    } else {
                        // cell is in halo
                        search_in_pair_of_cells(i,j,k,
                                                s1,s2,s3,
                                                grid1, grid2,
                                                bon,dist_vec,
                                                nlist.wrapped[i1] && is_periodic);
                        search_in_pair_of_cells(s1,s2,s3,
                                                i,j,k,
                                                grid1, grid2,
                                                bon,dist_vec,
                                                nlist.wrapped[i1] && is_periodic);
                    }


                }

            }
        }
    }
}

