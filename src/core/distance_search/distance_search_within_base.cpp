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



#include "distance_search_within_base.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;



void Distance_search_within_base::do_search()
{
    // Plan of search. Contain pairs of cells to search between
    // If two cells are the same then search inside the cell
    vector< Matrix<int,3,2> > plan;
    make_search_plan(plan);

    // Prepare for searching
    result->clear();

    // See if we need parallelization
    int nt = std::min(plan.size(), size_t(std::thread::hardware_concurrency()));

    if(nt==1){
        // Serial
        for(auto& pair: plan){
            search_planned_pair(pair.col(0),pair.col(1), *result);
        }

    } else {
        // Thread parallel
        vector<vector<int>> res_buf(nt);

        vector<thread> threads;
        int chunk = floor(plan.size()/nt);

        for(int i=0;i<nt;++i){
            int b = chunk*i;
            int e = (i<nt-1) ? chunk*(i+1) : plan.size();

            // Launch thread
            threads.emplace_back(
            [&,b,e,i](){
                for(int j=b; j<e; ++j){
                    search_planned_pair(plan[j].col(0), plan[j].col(1), res_buf[i]);
                }
            }
            );
        }

        // Wait for threads
        for(auto& t: threads) t.join();

        // Collect results
        for(int i=0;i<nt;++i){
            copy(begin(res_buf[i]),end(res_buf[i]),back_inserter(*result));
        }
    }
}

// grid1 is the source which is searched
void Distance_search_within_base::search_between_cells(Vector3i_const_ref c1,
                                                       Vector3i_const_ref c2,
                                                       const Grid &grid1,
                                                       const Grid &grid2,
                                                       std::vector<int> &res_buffer)
{
    const Grid_cell& cell1 = grid1.cell(c1);
    const Grid_cell& cell2 = grid2.cell(c2);

    // The cells could be not adjucent if one of them is wrapped around
    // but this only happens in peridic case and is treated automatically

    int N1 = cell1.size();
    int N2 = cell2.size();

    if(N1*N2==0) return; // Nothing to do

    float cutoff2 = cutoff*cutoff;

    if(is_periodic){

        for(int i1=0;i1<N1;++i1){
            Vector3f p = cell1.get_coord(i1); // Coord of point in grid1
            for(int i2=0;i2<N2;++i2){
                float d = box.distance_squared(cell2.get_coord(i2),p,periodic_dims);
                if(d<=cutoff2){
                    res_buffer.push_back(cell1.get_index(i1));
                }
            }
        }

    } else {

        for(int i1=0;i1<N1;++i1){
            Vector3f p = cell1.get_coord(i1); // Coord of point in grid1
            for(int i2=0;i2<N2;++i2){
                float d = (cell2.get_coord(i2)-p).squaredNorm();
                if(d<=cutoff2){                    
                    res_buffer.push_back(cell1.get_index(i1));
                }
            }
        }

    }
}

void Distance_search_within_base::search_planned_pair(Vector3i_const_ref c1, Vector3i_const_ref c2, std::vector<int> &res_buffer)
{
    search_between_cells(c1,c2,grid1,grid2,res_buffer);
    if(c1!=c2) search_between_cells(c1,c2,grid2,grid1,res_buffer);
}

