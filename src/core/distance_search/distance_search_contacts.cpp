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


#include "distance_search_contacts.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;

void Distance_search_contacts::do_search(){
    int i,j,k,i1;
    int nlist_size;

    // Search
    result_pairs->clear();
    if(result_distances) result_distances->clear();

    // Init visited cells array
    for(i=0;i<NgridX;++i)
        for(j=0;j<NgridY;++j)
            for(k=0;k<NgridZ;++k){
                visited[i][j][k] = false;
            }

    // See if we need parallelization
    int max_N, max_dim;
    Vector3i dims(NgridX,NgridY,NgridZ);
    max_N = dims.maxCoeff(&max_dim);

    int nt = std::min(max_N, int(std::thread::hardware_concurrency()));

    if(nt==1){
        deque<Vector2i> _bon;
        deque<float> _dist_vec;
        deque<float>* _dist_vec_ptr;
        _dist_vec_ptr = result_distances ? &_dist_vec : nullptr;

        do_part(0,0,NgridX,_bon,_dist_vec_ptr);

        // Collect results
        result_pairs->reserve(_bon.size());

        copy(_bon.begin(),_bon.end(),back_inserter(*result_pairs));

        if(result_distances){
            result_distances->reserve(_dist_vec.size());
            copy(_dist_vec.begin(),_dist_vec.end(),back_inserter(*result_distances));
        }

    } else {
        // Parallel searching

        // Determine parts for each thread
        vector<int> b(nt),e(nt);
        int cur=0;
        for(int i=0;i<nt-1;++i){
            b[i]=cur;
            cur += dims(max_dim)/nt;
            e[i]=cur;
        }
        b[nt-1]=cur;
        e[nt-1]=dims(max_dim);

        // Prepare arrays per each thread
        vector< deque<Vector2i> > _bon(nt);
        vector< deque<float> > _dist_vec(nt);
        vector< deque<float>* > _dist_vec_ptr(nt);
        for(int i=0;i<nt;++i) _dist_vec_ptr[i] = result_distances ? &_dist_vec[i] : nullptr;

        // Launch threads
        vector<thread> threads;
        for(int i=0;i<nt;++i){
            threads.emplace_back(
                                   std::bind(
                                       &Distance_search_contacts::do_part,
                                       this,
                                       max_dim,
                                       b[i],
                                       e[i],
                                       ref(_bon[i]),
                                       ref(_dist_vec_ptr[i])
                                   )                                
                             );
        }

        // Wait for threads
        for(auto& t: threads) t.join();

        // Collect results
        int sz = 0;
        for(int i=0;i<nt;++i) sz+= _bon[i].size();

        result_pairs->reserve(sz);
        for(int i=0;i<nt;++i){
            copy(_bon[i].begin(),_bon[i].end(),back_inserter(*result_pairs));
        }

        if(result_distances){
            result_distances->reserve(sz);
            for(int i=0;i<nt;++i)
                copy(_dist_vec[i].begin(),_dist_vec[i].end(),back_inserter(*result_distances));
        }
    }
}


void Distance_search_contacts::search_in_pair_of_cells(int x1, int y1, int z1, // cell 1
                             int x2, int y2, int z2, // cell 2
                             Grid& grid1,
                             Grid& grid2,
                             deque<Vector2i>& bon,
                             deque<float>* dist_vec, bool is_periodic)
{
    int N1,N2,ind1,ind2,i1,i2;
    float d;
    float cutoff2 = cutoff*cutoff;

    N1 = grid1.cell(x1,y1,z1).size();
    N2 = grid2.cell(x2,y2,z2).size();

    if(N1*N2==0) return; // Nothing to do

    const vector<Grid_element>& v1 = grid1.cell(x1,y1,z1);
    const vector<Grid_element>& v2 = grid2.cell(x2,y2,z2);

    if(is_periodic){

        for(i1=0;i1<N1;++i1){
            Vector3f* p = v1[i1].coor_ptr; // Coord of point in grid1
            for(i2=0;i2<N2;++i2){
                d = box.distance_squared(*(v2[i2].coor_ptr),*p);
                if(d<=cutoff2){
                    ind1 = v1[i1].index; //index
                    ind2 = v2[i2].index; //index
                    bon.emplace_back(ind1,ind2);
                    if(dist_vec) dist_vec->push_back(sqrt(d));
                }
            }
        }

    } else {

        for(i1=0;i1<N1;++i1){
            Vector3f* p = v1[i1].coor_ptr; // Coord of point in grid1
            for(i2=0;i2<N2;++i2){
                d = (*(v2[i2].coor_ptr)-*p).squaredNorm();
                if(d<=cutoff2){
                    ind1 = v1[i1].index; //index
                    ind2 = v2[i2].index; //index
                    bon.emplace_back(ind1,ind2);
                    if(dist_vec) dist_vec->push_back(sqrt(d));
                }
            }
        }

    }
}

