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



#include "distance_search_base.h"
#include "search_utils.h"
#include "pteros/core/pteros_error.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

// Non-periodic variant
void Distance_search_base::set_grid_size(const Vector3f &min, const Vector3f &max)
{
    Vector3f extents = max-min;    

    // Cell size should be >= cutoff for all dimentions
    NgridX = floor(extents(0)/cutoff);
    NgridY = floor(extents(1)/cutoff);
    NgridZ = floor(extents(2)/cutoff);

    if(NgridX<1) NgridX = 1;
    if(NgridY<1) NgridY = 1;
    if(NgridZ<1) NgridZ = 1;
}

// Periodic variant
void Distance_search_base::set_grid_size(const Periodic_box &box)
{
    Vector3f extents;
    extents(0) = box.box_to_lab(box.get_vector(0))(0);
    extents(1) = box.box_to_lab(box.get_vector(1))(1);
    extents(2) = box.box_to_lab(box.get_vector(2))(2);

    // Cell size should be >= cutoff for all dimentions
    NgridX = floor(extents(0)/cutoff);
    NgridY = floor(extents(1)/cutoff);
    NgridZ = floor(extents(2)/cutoff);

    if(NgridX<1) NgridX = 1;
    if(NgridY<1) NgridY = 1;
    if(NgridZ<1) NgridZ = 1;
}

void Distance_search_base::make_search_plan(vector<Matrix<int,3,2>>& plan){
    // Number of pairs to search:
    // X*Y*Z single cells
    // X*(Y*Z) pairs along X
    // Y*(X*Z) pairs along Y
    // Z*(X*Y) pairs along Z
    // = 4*X*Y*Z
    // In non-periodic variant less then this but we reserve memory for periodic variant anyway for simplicity
    plan.reserve(NgridX*NgridY*NgridZ*8);

    // Fill the plan
    // Go to the right of the current cell on each axis
    // In non-periodic variant skip wrapping
    Matrix<int,3,2> pair;
    for(int x=0;x<NgridX;++x){
        if(x==NgridX-1 && !periodic_dims(0)) continue;

        for(int y=0;y<NgridY;++y){
            if(y==NgridY-1 && !periodic_dims(1)) continue;

            for(int z=0;z<NgridZ;++z){
                if(z==NgridZ-1 && !periodic_dims(2)) continue;

                // We make a step forward by +1 in each direction
                for(int i1=0;i1<2;++i1){
                    for(int i2=0;i2<2;++i2){
                        for(int i3=0;i3<2;++i3){
                            pair << x,(x+i1) % NgridX,
                                    y,(y+i2) % NgridY,
                                    z,(z+i3) % NgridZ;
                            plan.push_back(pair);
                        } //i3
                    } //i2
                } //i1

            } //z
        } //y
    } //x
}

void Distance_search_base::create_grid(const Selection &sel)
{
    if(!is_periodic){
        sel.minmax(min,max);
        set_grid_size(min,max);
    } else {
        // Check if we have periodicity
        if(!box.is_periodic())
            throw Pteros_error("Asked for pbc in distance search, but there is no periodic box!");
        // Set dimensions of the current unit cell
        set_grid_size(box);
    }

    // Allocate one grid
    grid1.resize(NgridX,NgridY,NgridZ);
}

void Distance_search_base::create_grids(const Selection &sel1, const Selection &sel2)
{
    if(!is_periodic){
        // Get the minmax of each selection
        Vector3f min1,min2,max1,max2;

        sel1.minmax(min1,max1);
        sel2.minmax(min2,max2);        

        // Find bounding box
        for(int i=0;i<3;++i){
            overlap_1d(min1(i),max1(i),min2(i),max2(i),min(i),max(i));
            // If no overlap just exit
            if(max(i)==min(i)) return;
        }

        // Add a "halo: of size cutoff to bounding box
        min.array() -= cutoff;
        max.array() += cutoff;

        set_grid_size(min,max);

    } else {
        // Check if we have periodicity
        if(!box.is_periodic())
            throw Pteros_error("Asked for pbc in distance search, but there is no periodic box!");
        // Set dimensions of the current unit cell
        set_grid_size(box);
    }

    // Allocate both grids
    grid1.resize(NgridX,NgridY,NgridZ);
    grid2.resize(NgridX,NgridY,NgridZ);
}
