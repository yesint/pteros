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
        coor = sel.xyz_ptr(i);

        n1 = floor(NX*((*coor)(0)-min(0))/(max(0)-min(0)));
        if(n1<0 || n1>=NX) continue;

        n2 = floor(NY*((*coor)(1)-min(1))/(max(1)-min(1)));
        if(n2<0 || n2>=NY) continue;

        n3 = floor(NZ*((*coor)(2)-min(2))/(max(2)-min(2)));
        if(n3<0 || n3>=NZ) continue;

        if(abs_index){
            cell(n1,n2,n3).emplace_back(sel.index(i),coor);
        } else {
            cell(n1,n2,n3).emplace_back(i,coor);
        }
    }
}

void Grid::populate_periodic(const Selection &sel, bool abs_index)
{
    populate_periodic(sel, sel.box(), abs_index);
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
        coor = sel.xyz(i);
        // See if atom i is in box and wrap if needed
        if( !box.in_box(coor) ){
            box.wrap_point(coor);
            wrapped_atoms.push_back(coor);
            ptr = &*wrapped_atoms.rbegin();
        } else {
            ptr = sel.xyz_ptr(i);
        }

        // Now we are sure that coor is wrapped
        // Get relative coordinates in box [0:1)
        coor = m_inv*coor;

        n1 = floor(NX*coor(0));
        n2 = floor(NY*coor(1));
        n3 = floor(NZ*coor(2));

        // if coor(i) is 1.000001 or -0.00001 due to numerucal errors correct manually
        if(n1>=NX)
            n1=NX-1;
        else if(n1<0)
            n1=0;

        if(n2>=NY)
            n2=NY-1;
        else if(n2<0)
            n2=0;

        if(n3>=NZ)
            n3=NZ-1;
        else if(n3<0)
            n3=0;

        // Assign to grid
        if(abs_index){
            cell(n1,n2,n3).emplace_back(sel.index(i),ptr);
        } else {
            cell(n1,n2,n3).emplace_back(i,ptr);
        }
    }
}


