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





#include "distance_search_within_base.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;

void DistanceSearchWithinBase::compute_chunk(int b, int e,
                                             std::unordered_set<int>& res_buf)
{
    PlannedPair pair;
    for(int ind=b;ind<e;++ind){
        // Get array position
        auto p = index_to_pos(ind);
        // Apply stencil
        for(int j=0;j<stencil.size();j+=2){
            pair.c1 = p + stencil[j];
            pair.c2 = p + stencil[j+1];
            if(process_neighbour_pair(pair)){
                search_planned_pair(pair, res_buf);
            }
        }
    }
}


void DistanceSearchWithinBase::do_search()
{
    // Prepare for searching
    result.clear();

    // See if we need parallelization
    int nt = std::min(size_t(Ngrid.prod()), size_t(std::thread::hardware_concurrency()));

    if(nt==1){
        // Serial
        compute_chunk(0,Ngrid.prod(),result);

    } else {
        // Thread parallel
        vector<unordered_set<int>> res_buf(nt);

        vector<thread> threads;
        int chunk = floor(Ngrid.prod()/nt);

        for(int i=0;i<nt;++i){
            int b = chunk*i;
            int e = (i<nt-1) ? chunk*(i+1) : Ngrid.prod();

            // Launch thread
            threads.emplace_back(
            [&,b,e,i](){
                compute_chunk(b,e, res_buf[i]);
            }
            );
        }

        // Wait for threads
        for(auto& t: threads) t.join();

        // Collect results
        for(int i=0;i<nt;++i){
            // Copy to result guarantied to be unique, but not sorted
            result.insert(begin(res_buf[i]),end(res_buf[i]));
        }
    }
}

// grid1 is the source which is searched
void DistanceSearchWithinBase::search_between_cells(Vector3i_const_ref c1,
                                                       Vector3i_const_ref c2,
                                                       Vector3i_const_ref wrapped,
                                                       std::unordered_set<int> &res_buffer)
{
    const GridCell& cell1 = grid1.cell(c1);
    const GridCell& cell2 = grid2.cell(c2);

    // The cells could be not adjucent if one of them is wrapped around
    // but this only happens in peridic case and is treated automatically

    int N1 = cell1.size();
    int N2 = cell2.size();

    if(N1*N2==0) return; // Nothing to do

    float cutoff2 = cutoff*cutoff;

    if((wrapped.array()>0).any()){

        for(int i1=0;i1<N1;++i1){
            Vector3f p = cell1.get_coord(i1); // Coord of point in grid1
            for(int i2=0;i2<N2;++i2){
                float d = box.distance_squared(cell2.get_coord(i2), p, wrapped);
                if(d<=cutoff2){
                    res_buffer.insert(cell1.get_index(i1));
                }
            }
        }

    } else {

        for(int i1=0;i1<N1;++i1){
            Vector3f p = cell1.get_coord(i1); // Coord of point in grid1
            for(int i2=0;i2<N2;++i2){                
                float d = (cell2.get_coord(i2)-p).squaredNorm();                
                if(d<=cutoff2){                    
                    res_buffer.insert(cell1.get_index(i1));
                }
            }
        }

    }
}

void DistanceSearchWithinBase::search_planned_pair(const PlannedPair &pair, std::unordered_set<int> &res_buffer)
{
    search_between_cells(pair.c1,pair.c2,pair.wrapped,res_buffer);
    if(pair.c1!=pair.c2){
        search_between_cells(pair.c2,pair.c1,pair.wrapped,res_buffer);
    }
}



