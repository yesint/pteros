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

using namespace std;
using namespace pteros;
using namespace Eigen;

void Distance_search_base::set_grid_size(const Vector3f &min, const Vector3f &max, int Natoms, const Periodic_box &box)
{

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
    //TODO: This need to be refactored to get rid of the while loop and to
    // compute optimal size in one operation.
    if(is_periodic) {

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

    } else { // No projection needed since there is no box

        while(dX < cutoff && NgridX>1){
            --NgridX;
            dX = (max(0)-min(0))/NgridX;
        }
        while(dY < cutoff && NgridY>1){
            --NgridY;
            dY = (max(1)-min(1))/NgridY;
        }
        while(dZ < cutoff && NgridZ>1){
            --NgridZ;
            dZ = (max(2)-min(2))/NgridZ;
        }

    }
}

void Distance_search_base::get_nlist(int i, int j, int k, Nlist_t &nlist)
{

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


