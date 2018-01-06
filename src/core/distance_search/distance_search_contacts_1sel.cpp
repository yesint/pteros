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


#include "distance_search_contacts_1sel.h"
#include "pteros/core/pteros_error.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;


Distance_search_contacts_1sel::Distance_search_contacts_1sel(float d, const Selection &sel,
                            std::vector<Eigen::Vector2i>& pairs,
                            bool absolute_index,
                            bool periodic,
                            std::vector<float>* dist_vec){
    result_pairs = &pairs;
    result_distances = dist_vec;
    cutoff = d;
    is_periodic = periodic;
    abs_index = absolute_index;
    box = sel.Box();

    create_grid(sel);
    if(is_periodic){
        grid1.populate_periodic(sel,box,abs_index);
    } else {
        grid1.populate(sel,min,max,abs_index);
    }
    do_search();
}

void Distance_search_contacts_1sel::create_grid(const Selection &sel)
{
    if(!is_periodic){
        // Get the minmax of selection
        sel.minmax(min,max);
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

    set_grid_size(min,max, sel.size(), box);
    // Allocate one grid
    grid1.resize(NgridX,NgridY,NgridZ);
    // Allocate visited array
    visited.resize( boost::extents[NgridX][NgridY][NgridZ] );
}

void Distance_search_contacts_1sel::do_part(int dim, int _b, int _e, std::deque<Vector2i> &bon, std::deque<float> *dist_vec)
{
    Vector3i b(0,0,0);
    Vector3i e(NgridX,NgridY,NgridZ);
    int dim_max = e(dim);
    b(dim)= _b;
    e(dim)= _e;
    int i,j,k,i1,nlist_size;
    Nlist_t nlist; // Local nlist

    float cutoff2 = cutoff*cutoff;

    for(i=b(0);i<e(0);++i){
        for(j=b(1);j<e(1);++j){
            for(k=b(2);k<e(2);++k){
                // Search in central cell
                //get_central_1(i,j,k, sel, bon, dist_vec);
                search_in_cell(i,j,k,bon,dist_vec,false);
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

                    // We only check for visited cells inside local part, not in the "halo"
                    if(    cell(dim)>=b(dim)
                        && cell(dim)<e(dim) ){
                        // cell is inside the partition
                        if( !visited[cell(0)][cell(1)][cell(2)] )
                            search_in_pair_of_cells(i,j,k,
                                                    cell(0),cell(1),cell(2),
                                                    grid1, grid1,
                                                    bon,dist_vec,
                                                    nlist.wrapped[i1] && is_periodic);
                    } else {
                        // cell is in halo
                        search_in_pair_of_cells(i,j,k,
                                                cell(0),cell(1),cell(2),
                                                grid1, grid1,
                                                bon,dist_vec,
                                                nlist.wrapped[i1] && is_periodic);
                    }


                }

            }
        }
    }
}

void Distance_search_contacts_1sel::search_in_cell(int x, int y, int z,
                    deque<Vector2i>& bon,
                    deque<float>* dist_vec,
                    bool is_periodic)
    {
    int N,ind1,ind2,i1,i2;
    float d;
    float cutoff2 = cutoff*cutoff;

    N = grid1.cell(x,y,z).size();

    if(N==0) return; // Nothing to do

    const vector<Grid_element>& v = grid1.cell(x,y,z);

    // Absolute or local index is filled during filling the grid before

    if(is_periodic){

        for(i1=0;i1<N-1;++i1){
            Vector3f* p = v[i1].coor_ptr; // Coord of point in grid1
            for(i2=i1+1;i2<N;++i2){
                d = box.distance_squared(*(v[i2].coor_ptr),*p);
                if(d<=cutoff2){
                    ind1 = v[i1].index; //index
                    ind2 = v[i2].index; //index
                    bon.push_back(Vector2i(ind1,ind2));
                    if(dist_vec) dist_vec->push_back(sqrt(d));
                }
            }
        }

    } else {

        for(i1=0;i1<N-1;++i1){
            Vector3f* p = v[i1].coor_ptr; // Coord of point in grid1
            for(i2=i1+1;i2<N;++i2){
                d = (*(v[i2].coor_ptr)-*p).squaredNorm();
                if(d<=cutoff2){
                    ind1 = v[i1].index; //index
                    ind2 = v[i2].index; //index
                    bon.push_back(Vector2i(ind1,ind2));
                    if(dist_vec) dist_vec->push_back(sqrt(d));
                }
            }
        }

    }
}

