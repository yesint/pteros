/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/

#include "pteros/core/grid_search.h"
#include "pteros/core/selection.h"
#include "pteros/core/pteros_error.h"
#include <vector>
#include <map>
#include <algorithm>
#include "time.h"
#include <iostream>
#include <thread>
#include <array>

using namespace std;
using namespace pteros;
using namespace Eigen;

// Get intersection of two 1d bars
void overlap_1d(float a1, float a2, float b1, float b2, float& res1, float& res2){
    res1 = res2 = 0.0;
    if(a1<b1){
        if(a2<b1){
            return; // No overlap
         } else { //a2>b1
            res1 = b1;
            if(a2<b2) res2=a2; else res2=b2;
        }
    } else { //a1>b1
        if(a1>b2){
            return; //No overlap
        } else { //a1<b2
            res1 = a1;
            if(a2>b2) res2=b2; else res2=a2;
        }
    }
}

Grid_searcher::Grid_searcher(float d, const Selection &sel,
                            std::vector<Eigen::Vector2i>& bon,
                            bool absolute_index,
                            bool periodic,
                            std::vector<float>* dist_vec){
    cutoff = d;
    is_periodic = periodic;    
    abs_index = absolute_index;
    box = sel.get_system()->Box(sel.get_frame());    

    create_grid(grid1,sel);
    if(is_periodic){
        grid1.populate_periodic(sel,box,abs_index);
    } else {
        grid1.populate(sel,min,max,abs_index);
    }
    do_search1(bon,dist_vec);
}

Grid_searcher::Grid_searcher(float d, const Selection &sel1, const Selection &sel2,
                            std::vector<Eigen::Vector2i>& bon,
                            bool absolute_index,
                            bool periodic,
                            std::vector<float>* dist_vec){
    cutoff = d;
    is_periodic = periodic;    
    abs_index = absolute_index;
    box = sel1.get_system()->Box(sel1.get_frame());    

    create_grid2(sel1,sel2);

    if(is_periodic){
        grid1.populate_periodic(sel1,box,abs_index);
        grid2.populate_periodic(sel2,box,abs_index);
    } else {
        grid1.populate(sel1,min,max,abs_index);
        grid2.populate(sel2,min,max,abs_index);
    }

    do_search2(bon,dist_vec);
}

Grid_searcher::Grid_searcher(){
}


void Grid_searcher::assign_to_grid(float d, const Selection &sel,
                    bool absolute_index,
                    bool periodic){
    cutoff = d;
    is_periodic = periodic;
    abs_index = absolute_index;
    box = sel.get_system()->Box(sel.get_frame());    

    create_grid(grid1,sel);

    if(is_periodic){
        grid1.populate_periodic(sel,box,abs_index);
    } else {
        grid1.populate(sel,min,max,abs_index);
    }

    p_sel = const_cast<Selection*>(&sel);
}

void Grid_searcher::do_search_within(vector<int>& bon, const Selection& src){
    // Array of atomic bools for used source points
    std::vector<atomic_wrapper<bool>> used(src.size());
    for(int i=0;i<used.size();++i) used[i].store(false);

    //------------
    // Search part
    //------------
    bon.clear();

    Vector3f coor1;

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

        //for(int i=0;i<nt;++i) cout << b[i] << ":" << e[i] << " ";
        //cout << endl;

        // Launch threads
        vector<thread> threads;
        for(int i=0;i<nt;++i){
            threads.push_back( thread(
                                   std::bind(
                                       &Grid_searcher::do_part_within,
                                       //&Grid_searcher::do_part_within,
                                       this,
                                       max_dim,
                                       b[i],
                                       e[i],
                                       ref(used)
                                   )
                                )
                             );
        }

        // Wait for threads
        for(auto& t: threads) t.join();

    } else {
        // Serial search, no need to launch separate thread
        do_part_within(max_dim,0,dims(max_dim),used);
    }

    // Convert used array to indexes
    if(abs_index){
        for(int i=0;i<used.size();++i)
            if(used[i].load()) bon.push_back(src.Index(i));
    } else {
        for(int i=0;i<used.size();++i)
            if(used[i].load()) bon.push_back(i);
    }
}

void Grid_searcher::search_within(Vector3f_const_ref coord, vector<int> &bon){
    // grid1 has parent selection, which is "src". Our point is "target"
    System tmp; // tmp system with just one atom with coordinates coord
    vector<Vector3f> crd{coord};
    vector<Atom> atm(1);
    tmp.atoms_add(atm,crd);
    auto target = tmp.select_all();
    // Allocate second grid of the same size
    grid2.resize(NgridX,NgridY,NgridZ);

    if(is_periodic){
        grid2.populate_periodic(target,box,abs_index);
    } else {
        grid2.populate(target,min,max,abs_index);
    }

    // Now search
    do_search_within(bon,*p_sel);
}


void Grid_searcher::search_within(const Selection &target, std::vector<int> &bon, bool include_self){
    // Allocate second grid of the same size
    grid2.resize(NgridX,NgridY,NgridZ);

    Selection* ptr;
    Selection noself; // Only used for noself variant

    if(!include_self){
        noself = *p_sel;
        noself.remove(target);
        ptr = &noself;
    } else {
        ptr = p_sel;
    }

    if(is_periodic){
        grid2.populate_periodic(target,box,abs_index);
    } else {
        grid2.populate(target,min,max,abs_index);
    }

    do_search_within(bon,*ptr);
}

void search_in_cell(int x, int y, int z,
                    Grid& grid,
                    vector<Vector2i>& bon,
                    vector<float>* dist_vec,
                    const Periodic_box& box,
                    float cutoff2,
                    bool is_periodic, bool global_index)
    {
    int N,ind1,ind2,i1,i2;
    float d;

    N = grid.cell(x,y,z).size();

    if(N==0) return; // Nothing to do

    const vector<Grid_element>& v = grid.cell(x,y,z);

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

void search_in_pair_of_cells(int x1, int y1, int z1, // cell 1
                             int x2, int y2, int z2, // cell 2
                             Grid& grid1,
                             Grid& grid2,
                             vector<Vector2i>& bon,
                             vector<float>* dist_vec,
                             const Periodic_box& box,
                             float cutoff2,
                             bool is_periodic, bool global_index)
{
    int N1,N2,ind1,ind2,i1,i2;
    float d;

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
                    bon.push_back(Vector2i(ind1,ind2));
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
                    bon.push_back(Vector2i(ind1,ind2));
                    if(dist_vec) dist_vec->push_back(sqrt(d));
                }
            }
        }

    }
}

void search_in_pair_of_cells_for_within(int sx, int sy, int sz, // src cell
                             int tx, int ty, int tz, // target cell                             
                             Grid& grid1,
                             Grid& grid2,
                             std::vector<atomic_wrapper<bool>>& used,
                             const Periodic_box& box,
                             float cutoff2, bool is_periodic)
{
    int Ns,Nt,ind,s,t;    
    float d;

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

void Grid_searcher::do_part_within(int dim, int _b, int _e,
                             std::vector<atomic_wrapper<bool> > &used
                             ){
    Vector3i b(0,0,0);
    Vector3i e(NgridX,NgridY,NgridZ);
    int dim_max = e(dim);
    b(dim)= _b;
    e(dim)= _e;

    int i,j,k,c,t,ind;
    Nlist_t nlist; // Local nlist
    nlist.data.reserve(27);
    nlist.wrapped.reserve(27);

    float cutoff2 = cutoff*cutoff;

    for(i=b(0);i<e(0);++i){
        for(j=b(1);j<e(1);++j){
            for(k=b(2);k<e(2);++k){

                // Search in central cell
                search_in_pair_of_cells_for_within(i,j,k, //src cell
                                        i,j,k, //target cell
                                        grid1, grid2,
                                        used, box, cutoff2, false);
                // Get nlist                
                get_nlist(i,j,k,nlist);

                // Cycle over nlist
                for(c=0;c<nlist.data.size();++c){
                    const Vector3i& cell = nlist.data[c];

                    search_in_pair_of_cells_for_within(i,j,k, //src cell
                                            cell(0),cell(1),cell(2), //target cell
                                            grid1, grid2,
                                            used, box, cutoff2,
                                            nlist.wrapped[c] && is_periodic);
                }

            }
        }
    }

}

// Search is around target, atoms from src are returned
Grid_searcher::Grid_searcher(float d,
                            const Selection &src,
                            const Selection &target,
                            std::vector<int>& bon,
                            bool include_self,
                            bool periodic){

    cutoff = d;
    is_periodic = periodic;
    abs_index = true; // Force absolute indexes!

    // Get current box
    box = src.get_system()->Box(src.get_frame());

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

    // Fill grids (forced absolute indexes above!)
    if(is_periodic){
        grid1.populate_periodic(src,box,abs_index);
        grid2.populate_periodic(target,box,abs_index);
    } else {
        grid1.populate(src,min,max,abs_index);
        grid2.populate(target,min,max,abs_index);
    }

    const Selection* ptr;
    Selection noself; // Only used for noself variant

    if(!include_self){
        noself = src;
        noself.remove(target);
        ptr = &noself;
    } else {
        ptr = &src;
    }

    do_search_within(bon,*ptr);
}


void Grid_searcher::set_grid_size(const Vector3f& min, const Vector3f& max,
                                  int Natoms, const Periodic_box& box){

    /*  Our grids should satisfy these equations:
        NgridX * NgridY * NgridZ = Natoms
        NgridX/NgridY = a/b
        NgridY/NgridZ = b/c
        NgridX/NgridZ = a/c
        This lead to the following:
    */

    NgridX = floor(std::pow(double(Natoms*(max(0)-min(0))*(max(0)-min(0))/
                ((max(1)-min(1))*(max(2)-min(2)))), double(1.0/3.0))) ;
    NgridY = floor(std::pow(double(Natoms*(max(1)-min(1))*(max(1)-min(1))/
                ((max(0)-min(0))*(max(2)-min(2)))), double(1.0/3.0))) ;
    NgridZ = floor(std::pow(double(Natoms*(max(2)-min(2))*(max(2)-min(2))/
                ((max(0)-min(0))*(max(1)-min(1)))), double(1.0/3.0))) ;

    if(NgridX==0) NgridX = 1;
    if(NgridY==0) NgridY = 1;
    if(NgridZ==0) NgridZ = 1;    

    // Real grid vectors:
    float dX = (max(0)-min(0))/NgridX;
    float dY = (max(1)-min(1))/NgridY;
    float dZ = (max(2)-min(2))/NgridZ;

    // See if some of lab extents smaller than cutoff
    /*
    if(dX<cutoff) NgridX = floor(extX/cutoff);
    if(dY<cutoff) NgridY = floor(extY/cutoff);
    if(dZ<cutoff) NgridZ = floor(extZ/cutoff);
    */

    // See if some of grid vectors projected to lab axes smaller than cutoff
    while(box.box_to_lab(Vector3f(dX,0.0,0.0))(0) < cutoff && NgridX>1){
        --NgridX;
        dX = (max(0)-min(0))/NgridX;
    }
    while(box.box_to_lab(Vector3f(0.0,dY,0.0))(1) < cutoff && NgridY>1){
        --NgridY;
        dY = (max(1)-min(1))/NgridY;
    }
    while(box.box_to_lab(Vector3f(0.0,0.0,dZ))(2) < cutoff && NgridZ>1){
        --NgridZ;
        dZ = (max(2)-min(2))/NgridZ;
    }

    // See if some of lab extents larger than 2*cutoff
    /*
    if(dX>2.0*cutoff) NgridX = floor(extX/(2.0*cutoff));
    if(dY>2.0*cutoff) NgridY = floor(extY/(2.0*cutoff));
    if(dZ>2.0*cutoff) NgridZ = floor(extZ/(2.0*cutoff));
    */

    // See if some of grid vectors projected to lab axes larger than 2*cutoff
/*
    while(box.box_to_lab(Vector3f(dX,0.0,0.0))(0) > 2.0*cutoff){
        ++NgridX;
        dX = (max(0)-min(0))/NgridX;
    }
    while(box.box_to_lab(Vector3f(0.0,dY,0.0))(1) > 2.0*cutoff){
        ++NgridY;
        dY = (max(1)-min(1))/NgridY;
    }
    while(box.box_to_lab(Vector3f(0.0,0.0,dZ))(2) > 2.0*cutoff){
        ++NgridZ;
        dZ = (max(2)-min(2))/NgridZ;
    }
    */

}


void Grid_searcher::create_grid(Grid &grid, const Selection &sel)
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
    grid.resize(NgridX,NgridY,NgridZ);
    // Allocate visited array
    visited.resize( boost::extents[NgridX][NgridY][NgridZ] );
}


void Grid_searcher::create_grid2(const Selection &sel1, const Selection &sel2)
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


// Search over part of space. To be called in a thread.
void Grid_searcher::do_part1(int dim, int _b, int _e,
                             std::vector<Eigen::Vector2i>& bon,
                             std::vector<float>* dist_vec){

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
                search_in_cell(i,j,k, grid1,
                               bon,dist_vec,box,cutoff2,false,abs_index);
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
                                                    bon,dist_vec,box,cutoff2,
                                                    nlist.wrapped[i1] && is_periodic,
                                                    abs_index);
                    } else {
                        // cell is in halo                        
                        search_in_pair_of_cells(i,j,k,
                                                cell(0),cell(1),cell(2),
                                                grid1, grid1,
                                                bon,dist_vec,box,cutoff2,
                                                nlist.wrapped[i1] && is_periodic,
                                                abs_index);
                    }


                }

            }
        }
    }
}

// Search inside one selection
void Grid_searcher::do_search1(std::vector<Eigen::Vector2i>& bon,
                              std::vector<float>* dist_vec){
    int i,j,k,i1;
    int nlist_size;

    // Search
    bon.clear();
    if(dist_vec) dist_vec->clear();

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

    //int nt = std::min(max_N, int(std::thread::hardware_concurrency()));
    int nt=1;

    if(nt==1){
    //if(nt>0){
        // Serial searching
        do_part1(0,0,NgridX,bon,dist_vec);
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

        //for(int i=0;i<nt;++i) cout << b[i] << ":" << e[i] << " ";
        //cout << endl;

        // Prepare arrays per each thread
        vector< vector<Vector2i> > _bon(nt);
        vector< vector<float> > _dist_vec(nt);
        vector< vector<float>* > _dist_vec_ptr(nt);
        for(int i=0;i<nt;++i) _dist_vec_ptr[i] = dist_vec ? &_dist_vec[i] : nullptr;

        // Launch threads
        vector<thread> threads;
        for(int i=0;i<nt;++i){
            threads.push_back( thread(
                                   std::bind(
                                       &Grid_searcher::do_part1,
                                       this,
                                       max_dim,
                                       b[i],
                                       e[i],                                       
                                       ref(_bon[i]),
                                       ref(_dist_vec_ptr[i])
                                   )
                                )
                             );
        }

        // Wait for threads
        for(auto& t: threads) t.join();

        // Collect results
        for(int i=0;i<nt;++i){
            copy(_bon[i].begin(),_bon[i].end(), back_inserter(bon));
        }
        if(dist_vec){
            for(int i=0;i<nt;++i)
                copy(_dist_vec[i].begin(),_dist_vec[i].end(), back_inserter(*dist_vec));
        }
    }
}

// Search over part of space for two selections. To be called in a thread.
void Grid_searcher::do_part2(int dim, int _b, int _e,                             
                             std::vector<Eigen::Vector2i>& bon,
                             std::vector<float>* dist_vec){

    Vector3i b(0,0,0);
    Vector3i e(NgridX,NgridY,NgridZ);
    int dim_max = e(dim);
    b(dim)= _b;
    e(dim)= _e;
    int i,j,k,i1,nlist_size;
    Nlist_t nlist; // Local nlist

    float cutoff2 = cutoff*cutoff;

    int s1,s2,s3;

    for(i=b(0);i<e(0);++i){
        for(j=b(1);j<e(1);++j){
            for(k=b(2);k<e(2);++k){
                // Search in central cell
                // Central cell is always non-periodic
                search_in_pair_of_cells(i,j,k, i,j,k,
                                        grid1,grid2,
                                        bon,dist_vec,box,cutoff2,
                                        false,abs_index);
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
                                                    grid1,grid2,
                                                    bon,dist_vec,box,cutoff2,
                                                    nlist.wrapped[i1] && is_periodic,
                                                    abs_index);
                            search_in_pair_of_cells(s1,s2,s3,
                                                    i,j,k,
                                                    grid1,grid2,
                                                    bon,dist_vec,box,cutoff2,
                                                    nlist.wrapped[i1] && is_periodic,
                                                    abs_index);
                        }
                    } else {
                        // cell is in halo                        
                        search_in_pair_of_cells(i,j,k,
                                                s1,s2,s3,
                                                grid1,grid2,
                                                bon,dist_vec,box,cutoff2,
                                                nlist.wrapped[i1] && is_periodic,
                                                abs_index);
                        search_in_pair_of_cells(s1,s2,s3,
                                                i,j,k,
                                                grid1,grid2,
                                                bon,dist_vec,box,cutoff2,
                                                nlist.wrapped[i1] && is_periodic,
                                                abs_index);
                    }


                }

            }
        }
    }
}

// Search between two selections
void Grid_searcher::do_search2(std::vector<Eigen::Vector2i>& bon,
                              std::vector<float>* dist_vec){
    int i,j,k,nlist_size,i1;
    int n1,n2,n3;

    // Search
    bon.clear();    
    if(dist_vec) dist_vec->clear();

    // Init visited cells array
    for(i=0;i<NgridX;++i)
        for(j=0;j<NgridY;++j)
            for(k=0;k<NgridZ;++k)
                visited[i][j][k] = false;

    // See if we need parallelization
    int max_N, max_dim;
    Vector3i dims(NgridX,NgridY,NgridZ);
    max_N = dims.maxCoeff(&max_dim);

    int nt = std::min(max_N, int(std::thread::hardware_concurrency()));

    if(nt==1){
    //if(nt>0){
        // Serial search
        do_part2(0,0,NgridX,bon,dist_vec);
    } else {
        // Search in parallel
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

        //for(int i=0;i<nt;++i) cout << b[i] << ":" << e[i] << " ";
        //cout << endl;

        // Prepare arrays per each thread
        vector< vector<Vector2i> > _bon(nt);
        vector< vector<float> > _dist_vec(nt);
        vector< vector<float>* > _dist_vec_ptr(nt);
        for(int i=0;i<nt;++i) _dist_vec_ptr[i] = dist_vec ? &_dist_vec[i] : nullptr;

        // Launch threads
        vector<thread> threads;
        for(int i=0;i<nt;++i){
            threads.push_back( thread(
                                   std::bind(
                                       &Grid_searcher::do_part2,
                                       this,
                                       max_dim,
                                       b[i],
                                       e[i],
                                       ref(_bon[i]),
                                       ref(_dist_vec_ptr[i])
                                   )
                                )
                             );
        }

        // Wait for threads
        for(auto& t: threads) t.join();

        // Collect results
        for(int i=0;i<nt;++i){
            copy(_bon[i].begin(),_bon[i].end(), back_inserter(bon));
        }
        if(dist_vec){
            for(int i=0;i<nt;++i)
                copy(_dist_vec[i].begin(),_dist_vec[i].end(), back_inserter(*dist_vec));
        }
    }
}

void Grid_searcher::get_nlist(int i,int j,int k, Nlist_t& nlist){

    nlist.clear();

    Vector3i coor;

    if(!is_periodic){
        int c1,c2,c3;
        // Non-periodic variant
        for(c1=-1; c1<=1; ++c1){
            coor(0) = i+c1;
            if(coor(0)<0 || coor(0)>=NgridX) continue; // Bounds check
            for(c2=-1; c2<=1; ++c2){
                coor(1) = j+c2;
                if(coor(1)<0 || coor(1)>=NgridY) continue; // Bounds check
                for(c3=-1; c3<=1; ++c3){
                    coor(2) = k+c3;
                    if(coor(2)<0 || coor(2)>=NgridZ) continue; // Bounds check
                    //Exclude central cell
                    if(coor(0) == i && coor(1) == j && coor(2) == k ) continue;
                    // Add cell
                    nlist.append(coor);
                }
            }
        }
    } else {
        // Periodic variant
        int bX = 0, eX = 0;
        int bY = 0, eY = 0;
        int bZ = 0, eZ = 0;

        // If the number of cells in dimension is 2 this is a special case
        // when only one neighbour is need. Otherwise add both.
        if(NgridX>1) bX = -1;
        if(NgridY>1) bY = -1;
        if(NgridZ>1) bZ = -1;

        if(NgridX>2) eX = 1;
        if(NgridY>2) eY = 1;
        if(NgridZ>2) eZ = 1;

        int c1,c2,c3;        
        bool wrap1,wrap2,wrap3;

        for(c1 = bX; c1<=eX; ++c1){            
            wrap1 = false;
            coor(0) = i+c1;
            if(coor(0)==NgridX){ coor(0) = 0; wrap1=true; }
            if(coor(0)==-1){ coor(0) = NgridX-1; wrap1=true; }
            for(c2 = bY; c2<=eY; ++c2){
                wrap2 = false;
                coor(1) = j+c2;
                if(coor(1)==NgridY){ coor(1) = 0; wrap2=true; }
                if(coor(1)==-1){ coor(1) = NgridY-1; wrap2=true; }
                for(c3 = bZ; c3<=eZ; ++c3){
                    wrap3 = false;
                    coor(2) = k+c3;
                    if(coor(2)==NgridZ){ coor(2) = 0; wrap3=true; }
                    if(coor(2)==-1){ coor(2) = NgridZ-1; wrap3=true; }
                    //Exclude central cell
                    if(coor(0) == i && coor(1) == j && coor(2) == k) continue;
                    // Add cell
                    nlist.append(coor, wrap1||wrap2||wrap3);
                }
            }
        }
    }
}

void Nlist_t::clear()
{
    data.clear();
    wrapped.clear();
}

void Nlist_t::append(Vector3i_const_ref coor, bool wrap)
{
    data.push_back(coor);
    wrapped.push_back(wrap);
}
