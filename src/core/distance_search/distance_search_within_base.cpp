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


#include "distance_search_within_base.h"
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;

void Distance_search_within_base::used_to_result(vector<int>& res, bool include_self,
                                                 const Selection& src,
                                                 const Selection& target){
    res.clear();
    // Convert used to result
    if(include_self){
        if(abs_index){
            for(int i=0;i<used.size();++i)
                if(used[i].load()) res.push_back(src.Index(i));
        } else {
            for(int i=0;i<used.size();++i)
                if(used[i].load()) res.push_back(i);
        }
    } else {
        vector<int> dum;

        if(abs_index){
            for(int i=0;i<used.size();++i)
                if(used[i].load()) dum.push_back(src.Index(i));
        } else {
            for(int i=0;i<used.size();++i)
                if(used[i].load()) dum.push_back(i);
        }

        set_difference(dum.begin(),dum.end(),
                       target.index_begin(),target.index_end(),
                       back_inserter(res));
    }
}

void Distance_search_within_base::do_search(int sel_size)
{
    used.resize(sel_size);
    for(int i=0;i<used.size();++i) used[i].store(false);

    //------------
    // Search part
    //------------

    // See if we need parallelization
    int max_N, max_dim;
    Vector3i dims(NgridX,NgridY,NgridZ);
    max_N = dims.maxCoeff(&max_dim);

    int nt = std::min(max_N, int(std::thread::hardware_concurrency()));

    if(nt>1){
        // Parallel search

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

        // Launch threads
        vector<thread> threads;
        for(int i=0;i<nt;++i){
            threads.push_back( thread(
                                   std::bind(
                                       &Distance_search_within_base::do_part,
                                       this,
                                       max_dim,
                                       b[i],
                                       e[i]
                                   )
                                )
                             );
        }

        // Wait for threads
        for(auto& t: threads) t.join();

    } else {
        // Serial search, no need to launch separate thread
        do_part(max_dim,0,dims(max_dim));
    }

}


void Distance_search_within_base::do_part(int dim, int _b, int _e){
    Vector3i b(0,0,0);
    Vector3i e(NgridX,NgridY,NgridZ);
    int dim_max = e(dim);
    b(dim)= _b;
    e(dim)= _e;

    int i,j,k,c,t,ind;
    Nlist_t nlist; // Local nlist
    nlist.data.reserve(27);
    nlist.wrapped.reserve(27);

    for(i=b(0);i<e(0);++i){
        for(j=b(1);j<e(1);++j){
            for(k=b(2);k<e(2);++k){

                // Search in central cell
                search_in_pair_of_cells(i,j,k, //src cell
                                        i,j,k, //target cell
                                        false);
                // Get nlist
                get_nlist(i,j,k,nlist);

                // Cycle over nlist
                for(c=0;c<nlist.data.size();++c){
                    const Vector3i& cell = nlist.data[c];

                    search_in_pair_of_cells(i,j,k, //src cell
                                            cell(0),cell(1),cell(2), //target cell
                                            nlist.wrapped[c] && is_periodic);
                }

            }
        }
    }

}

void Distance_search_within_base::search_in_pair_of_cells(int sx, int sy, int sz, // src cell
                             int tx, int ty, int tz, // target cell
                             bool is_periodic)
{
    int Ns,Nt,ind,s,t;
    float d;
    float cutoff2 = cutoff*cutoff;

    Ns = grid1.cell(sx,sy,sz).size(); //src
    Nt = grid2.cell(tx,ty,tz).size(); //target

    if(Ns*Nt==0) return; // Nothing to do

    const vector<Grid_element>& sv = grid1.cell(sx,sy,sz);
    const vector<Grid_element>& tv = grid2.cell(tx,ty,tz);

    for(s=0;s<Ns;++s){
        ind = sv[s].index; // Local index here
        // Skip already used source points
        if(used[ind].load()) continue;

        Vector3f* p = sv[s].coor_ptr; // Coord of source point

        if(is_periodic){
            for(t=0;t<Nt;++t){
                d = box.distance_squared(*(tv[t].coor_ptr),*p);
                if(d<=cutoff2){
                    used[ind].store(true);
                    break;
                }
            }
        } else {
            for(t=0;t<Nt;++t){
                d = (*(tv[t].coor_ptr)-*p).squaredNorm();
                if(d<=cutoff2){
                    used[ind].store(true);
                    break;
                }
            }
        }

    }
}

