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


#include "distance_search_contacts.h"
#include <thread>
#include <algorithm>
#include <execution>

using namespace std;
using namespace pteros;
using namespace Eigen;


void Distance_search_contacts::do_search()
{
    // Plan of search. Contain pairs of cells to search between
    // If two cells are the same then search inside the cell
    vector< Matrix<int,3,2> > plan;
    make_search_plan(plan);

    // Prepare for searching
    pairs->clear();
    distances->clear();

    // See if we need parallelization    
    int nt = std::min(plan.size(), size_t(std::thread::hardware_concurrency()));

    if(nt==1){
        // Serial
        for(auto& pair: plan){
            search_planned_pair(pair.col(0),pair.col(1), *pairs, *distances);
        }

    } else {
        // Thread parallel
        vector<vector<Vector2i>> pairs_buf(nt);
        vector<vector<float>>    dist_buf(nt);

        vector<thread> threads;
        int chunk = floor(plan.size()/nt);

        for(int i=0;i<nt;++i){
            int b = chunk*i;
            int e = (i<nt-1) ? chunk*(i+1) : plan.size();

            // Launch thread
            threads.emplace_back(
            [&,b,e,i](){
                for(int j=b; j<e; ++j){
                    search_planned_pair(plan[j].col(0), plan[j].col(1), pairs_buf[i], dist_buf[i]);
                }
            }
            );
        }

        // Wait for threads
        for(auto& t: threads) t.join();

        // Collect results
        for(int i=0;i<nt;++i){
            copy(begin(pairs_buf[i]),end(pairs_buf[i]),back_inserter(*pairs));
            copy(begin(dist_buf[i]),end(dist_buf[i]),back_inserter(*distances));
        }
    }
}


void Distance_search_contacts::search_between_cells(Vector3i_const_ref c1,
                                                Vector3i_const_ref c2,
                                                const Grid& grid1,
                                                const Grid& grid2,
                                                std::vector<Eigen::Vector2i>& pairs_buffer,
                                                std::vector<float>& distances_buffer)
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
                    pairs_buffer.emplace_back(cell1.get_index(i1),
                                              cell2.get_index(i2));
                    distances_buffer.push_back(sqrt(d));
                }
            }
        }

    } else {

        for(int i1=0;i1<N1;++i1){
            Vector3f p = cell1.get_coord(i1); // Coord of point in grid1
            for(int i2=0;i2<N2;++i2){
                float d = (cell2.get_coord(i2)-p).squaredNorm();
                if(d<=cutoff2){
                    pairs_buffer.emplace_back(cell1.get_index(i1),
                                              cell2.get_index(i2));
                    distances_buffer.push_back(sqrt(d));
                }
            }
        }

    }
}

void Distance_search_contacts::search_inside_cell(Vector3i_const_ref c,
                                                  const Grid &grid,
                                                  std::vector<Vector2i> &pairs_buffer,
                                                  std::vector<float> &distances_buffer)
{
    const Grid_cell& cell = grid.cell(c);

    int N = cell.size();

    if(N==0) return; // Nothing to do

    float cutoff2 = cutoff*cutoff;

    if(is_periodic){

        for(int i1=0;i1<N-1;++i1){
            Vector3f p = cell.get_coord(i1); // Coord of point in grid1
            for(int i2=i1+1;i2<N;++i2){
                float d = box.distance_squared(cell.get_coord(i2),p,periodic_dims);
                if(d<=cutoff2){
                    pairs_buffer.emplace_back(cell.get_index(i1),
                                              cell.get_index(i2));
                    distances_buffer.push_back(sqrt(d));
                }
            }
        }

    } else {

        for(int i1=0;i1<N-1;++i1){
            Vector3f p = cell.get_coord(i1); // Coord of point in grid1
            for(int i2=i1+1;i2<N;++i2){
                float d = (cell.get_coord(i2)-p).squaredNorm();
                if(d<=cutoff2){
                    pairs_buffer.emplace_back(cell.get_index(i1),
                                              cell.get_index(i2));
                    distances_buffer.push_back(sqrt(d));
                }
            }
        }

    }
}




