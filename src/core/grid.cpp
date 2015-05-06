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

#include "pteros/core/grid.h"
#include "pteros/core/pteros_error.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

void Grid::clear()
{    
    int i,j,k;
    for(i=0;i<data.shape()[0];++i)
        for(j=0;j<data.shape()[1];++j)
            for(k=0;k<data.shape()[2];++k){
                data[i][j][k].clear();
        }

    wrapped_atoms.clear();
}

void Grid::resize(int X, int Y, int Z)
{
    data.resize( boost::extents[X][Y][Z] );
    clear();
}

void Grid::populate(const Selection &sel, bool abs_index)
{
    Vector3f min,max;
    sel.minmax(min,max);
    populate(sel,min,max,abs_index);
}

void Grid::populate(const Selection &sel, Vector3f_const_ref min, Vector3f_const_ref max, bool abs_index)
{
    int Natoms = sel.size();
    int NX = data.shape()[0];
    int NY = data.shape()[1];
    int NZ = data.shape()[2];
    int n1,n2,n3;

    // Non-periodic variant
    Vector3f* coor;
    for(int i=0;i<Natoms;++i){
        // Get coordinates of atom
        coor = sel.XYZ_ptr(i);

        n1 = floor(NX*((*coor)(0)-min(0))/(max(0)-min(0)));
        if(n1<0 || n1>=NX) continue;

        n2 = floor(NY*((*coor)(1)-min(1))/(max(1)-min(1)));
        if(n2<0 || n2>=NY) continue;

        n3 = floor(NZ*((*coor)(2)-min(2))/(max(2)-min(2)));
        if(n3<0 || n3>=NZ) continue;

        if(abs_index){
            cell(n1,n2,n3).push_back(Grid_element(sel.Index(i),coor));
        } else {
            cell(n1,n2,n3).push_back(Grid_element(i,coor));
        }
    }
}

void Grid::populate_periodic(const Selection &sel, bool abs_index)
{
    populate_periodic(sel, sel.Box(), abs_index);
}

void Grid::populate_periodic(const Selection &sel, const Periodic_box &box, bool abs_index)
{
    int Natoms = sel.size();
    int NX = data.shape()[0];
    int NY = data.shape()[1];
    int NZ = data.shape()[2];
    int n1,n2,n3;

    // Periodic variant
    Vector3f coor;
    Vector3f* ptr;
    Matrix3f m_inv = box.get_inv_matrix();

    for(int i=0;i<Natoms;++i){
        coor = sel.XYZ(i);
        // See if atom i is in box and wrap if needed
        if( !box.in_box(coor) ){
            box.wrap_point(coor);
            wrapped_atoms.push_back(coor);
            ptr = &*wrapped_atoms.rbegin();
        } else {
            ptr = sel.XYZ_ptr(i);
        }

        // Now we are sure that coor is wrapped
        // Get relative coordinates in box [0:1)
        coor = m_inv*coor;

        n1 = floor(NX*coor(0));
        n2 = floor(NY*coor(1));
        n3 = floor(NZ*coor(2));

        // if coor(i) is 1.000001 or -0.00001 due to numerucal errors correct manually
        if(n1>=NX) n1=NX-1;
        if(n1<0) n1=0;
        if(n2>=NY) n2=NY-1;
        if(n2<0) n2=0;
        if(n3>=NZ) n3=NZ-1;
        if(n3<0) n3=0;

        // Assign to grid
        if(abs_index){
            cell(n1,n2,n3).push_back(Grid_element(sel.Index(i),ptr));
        } else {
            cell(n1,n2,n3).push_back(Grid_element(i,ptr));
        }
    }
}
