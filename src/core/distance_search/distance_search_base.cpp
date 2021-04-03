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


#include "distance_search_base.h"
#include "pteros/core/pteros_error.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

// Non-periodic variant
void DistanceSearchBase::set_grid_size(const Vector3f &min, const Vector3f &max)
{
    Vector3f extents = max-min;    

    // Cell size should be >= cutoff for all dimentions
    Ngrid(0) = floor(extents(0)/cutoff);
    Ngrid(1) = floor(extents(1)/cutoff);
    Ngrid(2) = floor(extents(2)/cutoff);

    if(Ngrid(0)<1) Ngrid(0) = 1;
    if(Ngrid(1)<1) Ngrid(1) = 1;
    if(Ngrid(2)<1) Ngrid(2) = 1;
}

// Periodic variant
void DistanceSearchBase::set_grid_size(const PeriodicBox &box)
{
    Vector3f extents;
    extents(0) = box.box_to_lab(box.get_vector(0))(0);
    extents(1) = box.box_to_lab(box.get_vector(1))(1);
    extents(2) = box.box_to_lab(box.get_vector(2))(2);

    // Cell size should be >= cutoff for all dimentions
    Ngrid(0) = floor(extents(0)/cutoff);
    Ngrid(1) = floor(extents(1)/cutoff);
    Ngrid(2) = floor(extents(2)/cutoff);

    if(Ngrid(0)<1) Ngrid(0) = 1;
    if(Ngrid(1)<1) Ngrid(1) = 1;
    if(Ngrid(2)<1) Ngrid(2) = 1;
}

bool DistanceSearchBase::process_neighbour_pair(PlannedPair& pair){
    pair.wrapped.fill(0);
    for(int dim=0;dim<3;++dim){
        if(pair.c1(dim)==Ngrid(dim)){ // point beyond the right edge
            if(periodic_dims(dim)){
                pair.c1(dim) = pair.c1(dim) % Ngrid(dim); // Wrap this dimension
                pair.wrapped(dim) = 1;
            } else {
                return false; // don't use this pair
            }
        }

        if(pair.c2(dim)==Ngrid(dim)){ // point beyond the right edge
            if(periodic_dims(dim)){
                pair.c2(dim) = pair.c2(dim) % Ngrid(dim); // Wrap this dimension
                pair.wrapped(dim) = 1;
            } else {
                return false; // don't use this pair
            }
        }

        // Corner case: if cutoff>=0.5*extent then all pairs along this extent should
        // be forced periodic.
        if(periodic_dims(dim) && Ngrid(dim)<=2) pair.wrapped(dim) = 1;
    }

    return true; // use this pair
}


Vector3i DistanceSearchBase::index_to_pos(int i){
    Vector3i pos;
    // i = z+Nz*y+Nz*Ny*x
    pos(2) = i % Ngrid(2);
    pos(1) = (i/Ngrid(2)) % Ngrid(1);
    pos(0) = (i/Ngrid(2)/Ngrid(1)) % Ngrid(0);
    return pos;
}

/*
void Distance_search_base::make_search_plan(vector<Planned_pair>& plan){
    plan.reserve(Ngrid.prod()*stencil.size());

    Planned_pair pair;

    for(int x=0;x<Ngrid(0);++x){
        for(int y=0;y<Ngrid(1);++y){
            for(int z=0;z<Ngrid(2);++z){

                Vector3i p(x,y,z);
                for(int i=0;i<stencil.size();i+=2){
                    pair.c1 = p + stencil[i];
                    pair.c2 = p + stencil[i+1];
                    if(process_neighbour_pair(pair)) plan.push_back(pair);
                }

            } //z
        } //y
    } //x

}
*/

void DistanceSearchBase::create_grid(const Selection &sel)
{
    if(!is_periodic){
        sel.minmax(min,max);
        // add small margin to tolerate numeric errors
        max.array() += 1e-5;
        min.array() -= 1e-5;
        set_grid_size(min,max);
    } else {
        // Check if we have periodicity
        if(!box.is_periodic())
            throw PterosError("Asked for pbc in distance search, but there is no periodic box!");
        // Set dimensions of the current unit cell
        set_grid_size(box);
    }

    // Allocate one grid
    grid1.resize(Ngrid(0),Ngrid(1),Ngrid(2));
}


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



void DistanceSearchBase::create_grids(const Selection &sel1, const Selection &sel2)
{
    if(!is_periodic){
        // Get the minmax of each selection
        Vector3f min1,min2,max1,max2;

        sel1.minmax(min1,max1);
        sel2.minmax(min2,max2);

        // Find bounding box (intersection is tested with cutoff halo)
        for(int i=0;i<3;++i){
            overlap_1d(min1(i)-cutoff,max1(i)+cutoff,min2(i)-cutoff,max2(i)+cutoff,min(i),max(i));
            // If no overlap just exit
            if(max(i)==min(i)) return;
        }

        // In case of intersection add small margin to tolerate numeric errors
        max.array() += 1e-5;
        min.array() -= 1e-5;

        set_grid_size(min,max);

    } else {
        // Check if we have periodicity
        if(!box.is_periodic())
            throw PterosError("Asked for pbc in distance search, but there is no periodic box!");
        // Set dimensions of the current unit cell
        set_grid_size(box);
    }

    // Allocate both grids
    grid1.resize(Ngrid(0),Ngrid(1),Ngrid(2));
    grid2.resize(Ngrid(0),Ngrid(1),Ngrid(2));
}


