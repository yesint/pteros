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
#include <vector>
#include <map>
#include <algorithm>
#include <boost/multi_array.hpp>
#include "time.h"
#include <iostream>

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

Grid_searcher::Grid_searcher(float d, Selection& sel,
                            std::vector<Eigen::Vector2i>& bon,
                            bool absolute_index,
                            bool periodic,
                            std::vector<float>* dist_vec){
    cutoff = d;
    is_periodic = periodic;    
    abs_index = absolute_index;

    create_grid(grid1,sel);
    populate_grid(grid1,sel);
    do_search(sel,bon,dist_vec);
}

Grid_searcher::Grid_searcher(float d, Selection& sel1, Selection& sel2,
                            std::vector<Eigen::Vector2i>& bon,
                            bool absolute_index,
                            bool periodic,
                            std::vector<float>* dist_vec){
    cutoff = d;
    is_periodic = periodic;    
    abs_index = absolute_index;

    create_grid2(sel1,sel2);
    populate_grid(grid1,sel1);
    populate_grid(grid2,sel2);
    do_search(sel1,sel2,bon,dist_vec);
}

Grid_searcher::Grid_searcher(){
}


void Grid_searcher::assign_to_grid(float d, Selection& sel,
                    bool absolute_index,
                    bool periodic){
    cutoff = d;
    is_periodic = periodic;
    abs_index = absolute_index;
    create_grid(grid1,sel);
    populate_grid(grid1,sel);
    // Remember pointer to this selection
    p_sel = &sel;
}

void Grid_searcher::create_custom_grid(int nX, int nY, int nZ){
    NgridX = nX;
    NgridY = nY;
    NgridZ = nZ;

    grid1.resize( boost::extents[NgridX][NgridY][NgridZ] );
}

void Grid_searcher::fill_custom_grid(Selection sel, bool absolute_index){

    Vector3f box_dim = sel.get_system()->Box(sel.get_frame()).colwise().norm();
    min.fill(0.0);
    max = box_dim;
    // Grid vectors:
    dX = (max(0)-min(0))/NgridX;
    dY = (max(1)-min(1))/NgridX;
    dZ = (max(2)-min(2))/NgridX;

    is_periodic = true;
    abs_index = absolute_index;
    is_triclinic = false;
    populate_grid(grid1,sel);
}

vector<int>& Grid_searcher::cell_of_custom_grid(int x, int y, int z){
    return grid1[x][y][z];
}

void Grid_searcher::search_within(Vector3f &coor, vector<int> &bon){
    // Determine cell, which given point belongs to
    int n1,n2,n3,i,m1,m2,m3;

    bon.clear();

    // Get coordinates in triclinic basis if needed
    if(is_periodic && is_triclinic) coor = inv_basis_matr*coor;

    // Assign to grid
    n1 = floor((NgridX-1)*(coor(0)-min(0))/(max(0)-min(0)));
    n2 = floor((NgridY-1)*(coor(1)-min(1))/(max(1)-min(1)));
    n3 = floor((NgridZ-1)*(coor(2)-min(2))/(max(2)-min(2)));

    if(is_periodic){
        // If periodic and extends over the grid dimensions wrap it
        while(n1>=NgridX || n1<0)
            n1>=0 ? n1 %= NgridX : n1 = NgridX + n1%NgridX;
        while(n2>=NgridY || n2<0)
            n2>=0 ? n2 %= NgridY : n2 = NgridY + n2%NgridY;
        while(n3>=NgridZ || n3<0)
            n3>=0 ? n3 %= NgridZ : n3 = NgridZ + n3%NgridZ;
    } else {
        // In non-periodic variant discard point if doesn't fit into bounding box
        if(n1<0 || n1>=NgridX || n2<0 || n2>=NgridY || n3<0 || n3>=NgridZ) return;
    }

    float d;

    // Get neighbour list
    get_nlist(n1,n2,n3);
    // Add central cell to the list
    nlist.push_back(Vector3i(n1,n2,n3));
    int nlist_size = nlist.size();
    // Searh in all cells
    for(i=0;i<nlist_size;++i){
        m1 = nlist[i](0);
        m2 = nlist[i](1);
        m3 = nlist[i](2);
        int n = grid1[m1][m2][m3].size();
        for(int c=0;c<n;++c){            
            if(!is_periodic)
                d = (p_sel->XYZ(grid1[m1][m2][m3][c]) - coor).norm();
            else
                d = periodic_distance(p_sel->XYZ(grid1[m1][m2][m3][c]),coor);

            if(d<=cutoff){
                if(abs_index){
                    bon.push_back( p_sel->Index(grid1[m1][m2][m3][c]) );
                } else {
                    bon.push_back( grid1[m1][m2][m3][c] );
                }                
            }            

        }
    }
}


void Grid_searcher::search_within(Selection &target, std::vector<int> &bon, bool include_self){
    bon.clear();  

    int i,j,k,c,n1,n2,nlist_size,m1,m2,m3,N1,N2;
    float d;
    Vector3f coor1,coor2;

    // Allocate grid2 and populate it from target
    grid2.resize( boost::extents[NgridX][NgridY][NgridZ] );
    populate_grid(grid2,target);

    //cout << "Grid dimensions: " << NgridX << " " << NgridY << " " << NgridZ << endl;

    // Cycle over all cells of grid2
    for(i=0;i<NgridX;++i){
        for(j=0;j<NgridY;++j){
            for(k=0;k<NgridZ;++k){
                // Get number of atoms in current grid2 cell
                N2 = grid2[i][j][k].size();
                // If no atoms than just skip this cell
                if(N2==0) continue;
                // Get neighbour list
                get_nlist(i,j,k);
                // Add central cell to the list                
                nlist.push_back(Vector3i(i,j,k));

                nlist_size = nlist.size();
                // Cycle over neighbouring cells
                //cout << endl;
                for(c=0;c<nlist_size;++c){

                    m1 = nlist[c](0);
                    m2 = nlist[c](1);
                    m3 = nlist[c](2);

                    // Get number of atoms in neighbour grid1 cell
                    N1 = grid1[m1][m2][m3].size();
                    // Skip empty pairs
                    if(N1==0) continue;

                    // Cycle over N2 and N1                    
                    for(n2=0;n2<N2;++n2){
                        coor1 = target.XYZ(grid2[i][j][k][n2]);

                        for(n1=0;n1<N1;++n1){
                            // Skip already used points
                            if(grid1[m1][m2][m3][n1]<0) continue;
                            coor2 = p_sel->XYZ(grid1[m1][m2][m3][n1]);

                            if(!is_periodic)
                                d = (coor2 - coor1).norm();
                            else
                                d = periodic_distance(coor2,coor1);

                            if(d<=cutoff){
                                if(abs_index){
                                    bon.push_back( p_sel->Index(grid1[m1][m2][m3][n1]) );
                                } else {
                                    bon.push_back( grid1[m1][m2][m3][n1] );
                                }
                                // Mark atom in grid1 as already added
                                grid1[m1][m2][m3][n1] = -grid1[m1][m2][m3][n1];                                
                            }
                        }
                    }

                }
            }
        }
    }

    // Restore grid1 for possible later searches
    for(i=0;i<NgridX;++i)
        for(j=0;j<NgridY;++j)
            for(k=0;k<NgridZ;++k)
                for(n1=0;n1<grid1[i][j][k].size();++n1)
                    grid1[i][j][k][n1] = abs(grid1[i][j][k][n1]);

    if(include_self){
        // Add all target atoms to result
        copy(target.index_begin(),target.index_end(),back_inserter(bon));
    }

    sort(bon.begin(),bon.end());
    // Remove duplicates
    vector<int>::iterator it = std::unique(bon.begin(), bon.end());
    // Get rid of the tail with garbage
    bon.resize( it - bon.begin() );

    if(!include_self){
        vector<int> dum = bon;
        bon.clear();
        set_difference(dum.begin(),dum.end(),target.index_begin(),target.index_end(),back_inserter(bon));
    }
}

// Search is around target, atoms from src are returned
Grid_searcher::Grid_searcher(float d,
                            Selection& src,
                            Selection& target,
                            std::vector<int>& bon,
                            bool include_self,
                            bool absolute_index,
                            bool periodic){

    cutoff = d;
    is_periodic = periodic;
    abs_index = absolute_index;

    //------------
    // Grid creation part
    //------------

    // Get the minmax of each selection
    Vector3f min1,min2,max1,max2;
    int i;

    src.minmax(min1,max1);
    target.minmax(min2,max2);
    // Add a "halo: of size cutoff for each of them
    min1.array() -= cutoff;
    max1.array() += cutoff;
    min2.array() -= cutoff;
    max2.array() += cutoff;

    int effective_num = src.size() + target.size();

    // Get current box
    Matrix3f box = src.get_system()->Box(src.get_frame());
    if(!is_periodic){
        // Find true bounding box
        for(i=0;i<3;++i){
            overlap_1d(min1(i),max1(i),min2(i),max2(i),min(i),max(i));
            // If no overlap just exit
            if(max(i)==min(i)) return;
        }

        /*
        int num1 = 0, num2 = 0;
        Vector3f coor1;

        for(i=0; i<src.size(); ++i){
            coor1 = src.XYZ(i);
            if(   min(0)<coor1(0) && min(1)<coor1(1) && min(2)<coor1(2)
               && max(0)>coor1(0) && max(1)>coor1(1) && max(2)>coor1(2) ) ++num1;
        }

        for(i=0; i<target.size(); ++i){
            coor1 = src.XYZ(i);
            if(   min(0)<coor1(0) && min(1)<coor1(1) && min(2)<coor1(2)
               && max(0)>coor1(0) && max(1)>coor1(1) && max(2)>coor1(2) ) ++num2;
        }

        // See which selection have more atoms
        effective_num = std::max(num1,num2);
        //cout << effective_num << " " << src.size() << " " << target.size() << endl;
        */
    } else {
        // Set dimensions of the current unit cell
        min.fill(0.0);
        max = box_dim = box.colwise().norm(); //Also sets box dimensions

        // Compute basis conversion matrix if needed
        if(is_triclinic){
            make_inv_matr(box);
        }
    }

    set_grid_size(min,max, src.size()+target.size());
    //set_grid_size(min,max, effective_num);

    // Allocate both grids
    grid1.resize( boost::extents[NgridX][NgridY][NgridZ] );
    grid2.resize( boost::extents[NgridX][NgridY][NgridZ] );
    // Fill grids
    populate_grid(grid1,src);
    populate_grid(grid2,target);

    //------------
    // Search part
    //------------
    bon.clear();

    int j,k,c,n1,n2,nlist_size,m1,m2,m3,N1,N2,ind;

    Vector3f coor1;

    // Cycle over all cells of grid2 (target)
    for(i=0;i<NgridX;++i){
        for(j=0;j<NgridY;++j){
            for(k=0;k<NgridZ;++k){
                // Get number of atoms in current target cell
                N2 = grid2[i][j][k].size();
                // If no atoms than just skip this cell
                if(N2==0) continue;                

                // Matrix of pre-computed coordinates for target points in this cell
                // Saves access to grid and computing XYZ in place. Makes a big
                // difference for large cutoffs!
                MatrixXf pre(3,N2);
                for(n2=0;n2<N2;++n2){ //over target atoms
                    pre.col(n2) = target.XYZ(grid2[i][j][k][n2]);
                }


                // Get neighbour list
                get_nlist(i,j,k);
                // Add central cell to the list
                nlist.push_back(Vector3i(i,j,k));

                nlist_size = nlist.size();
                // Cycle over neighbouring cells                
                for(c=0;c<nlist_size;++c){

                    m1 = nlist[c](0);
                    m2 = nlist[c](1);
                    m3 = nlist[c](2);


                /*
                // Cycle over all accessible neighbouring cells
                for(m1=i-1;m1<=i+1;++m1)
                for(m2=j-1;m2<=j+1;++m2)
                for(m3=k-1;m3<=k+1;++m3)
                {
                    //Bounds check
                    if(m1<0 || m2<0 || m3<0 || m1>=NgridX || m2>=NgridY || m3>=NgridZ) continue;
                */


                    // Get number of atoms in neighbour grid1 cell
                    N1 = grid1[m1][m2][m3].size();

                    // Skip empty pairs
                    if(N1==0) continue;

                    // Cycle over N2 and N1                    
                    for(n1=0;n1<N1;++n1){

                        ind = grid1[m1][m2][m3][n1];
                        // Skip already used source points
                        if(ind<0) continue;

                        coor1 = src.XYZ(ind);

                        for(n2=0;n2<N2;++n2){ //over target atoms

                            if(!is_periodic)
                                d = (pre.col(n2) - coor1).norm();
                            else
                                d = periodic_distance(pre.col(n2), coor1);

                            if(d<=cutoff){
                                if(abs_index){
                                    bon.push_back( src.Index(ind) );
                                } else {
                                    bon.push_back( ind );
                                }
                                // Mark atom in grid1 as already added
                                grid1[m1][m2][m3][n1] = -ind;
                                // And break from cycle over n2 since atom is added already
                                break;
                            }
                        }
                    }


                    //--
                }
            }
        }
    }

    /*
    // Restore grid1 for possible later searches
    for(i=0;i<NgridX;++i)
        for(j=0;j<NgridY;++j)
            for(k=0;k<NgridZ;++k)
                for(n1=0;n1<grid1[i][j][k].size();++n1)
                    grid1[i][j][k][n1] = abs(grid1[i][j][k][n1]);
    */

    if(include_self){
        // Add all target atoms to result
        copy(target.index_begin(),target.index_end(),back_inserter(bon));
    }

    sort(bon.begin(),bon.end());
    // Remove duplicates
    vector<int>::iterator it = std::unique(bon.begin(), bon.end());
    // Get rid of the tail with garbage
    bon.resize( it - bon.begin() );

    if(!include_self){
        vector<int> dum = bon;
        bon.clear();
        set_difference(dum.begin(),dum.end(),target.index_begin(),target.index_end(),back_inserter(bon));
    }
}


void Grid_searcher::set_grid_size(const Vector3f& min, const Vector3f& max, int Natoms){
    /*  Our grids should satisfy these equations:
        NgridX * NgridY * NgridZ = Natoms
        NgridX/NgridY = a/b
        NgridY/NgridZ = b/c
        NgridX/NgridZ = a/c
        This lead to the following:
    */


    NgridX = floor(pow(Natoms*(max(0)-min(0))*(max(0)-min(0))/
                ((max(1)-min(1))*(max(2)-min(2))), 1.0/3.0)) ;
    NgridY = floor(pow(Natoms*(max(1)-min(1))*(max(1)-min(1))/
                ((max(0)-min(0))*(max(2)-min(2))), 1.0/3.0)) ;
    NgridZ = floor(pow(Natoms*(max(2)-min(2))*(max(2)-min(2))/
                ((max(0)-min(0))*(max(1)-min(1))), 1.0/3.0)) ;

    //NgridX = NgridY = NgridZ = pow(Natoms,1.0/3.0);

    // Coeff 3 is chosen empirically

    if(NgridX==0) NgridX = 1;
    if(NgridY==0) NgridY = 1;
    if(NgridZ==0) NgridZ = 1;    

    // Grid vectors:
    dX = (max(0)-min(0))/NgridX;
    dY = (max(1)-min(1))/NgridX;
    dZ = (max(2)-min(2))/NgridX;


    // See if some of grid vectors smaller then cutoff
    if(dX<cutoff){
        NgridX = floor((max(0)-min(0))/cutoff);
        dX = (max(0)-min(0))/NgridX;
    }
    if(dY<cutoff){
        NgridY = floor((max(1)-min(1))/cutoff);
        dY = (max(1)-min(1))/NgridY;
    }
    if(dZ<cutoff){
        NgridZ = floor((max(2)-min(2))/cutoff);
        dZ = (max(2)-min(2))/NgridZ;
    }


    //cout << "]]]]]]]]] " << NgridX << " " << NgridY << " " << NgridZ << endl;
}

void Grid_searcher::make_inv_matr(const Matrix3f& box){
    inv_basis_matr.col(0) = box.col(0).normalized();
    inv_basis_matr.col(1) = box.col(1).normalized();
    inv_basis_matr.col(2) = box.col(2).normalized();
    inv_basis_matr = inv_basis_matr.inverse().eval();
}

void Grid_searcher::create_grid(Grid_t& grid,Selection& sel){
    if(!is_periodic){
        // Get the minmax of selection
        sel.minmax(min,max);
        // Add a "halo: of size cutoff
        min.array() -= cutoff;
        max.array() += cutoff;
    } else {
        // Get current box
        Matrix3f box = sel.get_system()->Box(sel.get_frame());
        // Set dimensions of the current unit cell
        min.fill(0.0);
        max = box_dim = box.colwise().norm(); //Also sets box dimensions
        // Also compute basis conversion matrix
        if(is_triclinic){
            make_inv_matr(box);
        }
    }

    set_grid_size(min,max, sel.size());
    // Allocate one grid
    grid.resize( boost::extents[NgridX][NgridY][NgridZ] );
    // Allocate visited array
    visited.resize( boost::extents[NgridX][NgridY][NgridZ] );
}


void Grid_searcher::create_grid2(Selection& sel1, Selection& sel2){
    // Get the minmax of each selection
    Vector3f min1,min2,max1,max2;
    int i;

    sel1.minmax(min1,max1);
    sel2.minmax(min2,max2);
    // Add a "halo: of size cutoff for each of them
    min1.array() -= cutoff;
    max1.array() += cutoff;
    min2.array() -= cutoff;
    max2.array() += cutoff;

    // Get current box
    Matrix3f box = sel1.get_system()->Box(sel1.get_frame());
    if(!is_periodic)
        // Find true bounding box
        for(i=0;i<3;++i){
            overlap_1d(min1(i),max1(i),min2(i),max2(i),min(i),max(i));
            // If no overlap just exit
            if(max(i)==min(i)) return;
        }
    else {
        // Set dimensions of the current unit cell
        min.fill(0.0);
        max = box_dim = box.colwise().norm(); //Also sets box dimensions

        // Compute basis conversion matrix if needed
        if(is_triclinic){
            make_inv_matr(box);
        }
    }

    set_grid_size(min,max, sel1.size()+sel2.size());
    // Allocate both grids
    grid1.resize( boost::extents[NgridX][NgridY][NgridZ] );
    grid2.resize( boost::extents[NgridX][NgridY][NgridZ] );
    // Allocate visited array
    visited.resize( boost::extents[NgridX][NgridY][NgridZ] );
}

// General function for populating given grid from given selection
void Grid_searcher::populate_grid(Grid_t& grid, Selection& sel){
    int Natoms = sel.size();
    int n1,n2,n3,i,j,k;

    // Clear grid
    for(i=0;i<NgridX;++i)
        for(j=0;j<NgridY;++j)
            for(k=0;k<NgridZ;++k)
                grid[i][j][k].clear();

    // Assigning atoms to grid
    Vector3f coor;
    if(!is_periodic)
        // Non-periodic variant
        for(i=0;i<Natoms;++i){
            // Get coordinates of atom
            coor = sel.XYZ(i);

            n1 = floor((NgridX-0)*(coor(0)-min(0))/(max(0)-min(0)));
            if(n1<0 || n1>=NgridX) continue;

            n2 = floor((NgridY-0)*(coor(1)-min(1))/(max(1)-min(1)));
            if(n2<0 || n2>=NgridY) continue;

            n3 = floor((NgridZ-0)*(coor(2)-min(2))/(max(2)-min(2)));
            if(n3<0 || n3>=NgridZ) continue;

            grid[n1][n2][n3].push_back(i);
        }
    else {
        // Periodic variant        
        for(i=0;i<Natoms;++i){
            // Get coordinates of atom
            coor = sel.XYZ(i);
            // Get coordinates in triclinic basis if needed
            if(is_triclinic) coor = inv_basis_matr*coor;
            // Assign to non-periodic grid first
            n1 = floor((NgridX-0)*(coor(0)-min(0))/(max(0)-min(0)));
            n2 = floor((NgridY-0)*(coor(1)-min(1))/(max(1)-min(1)));
            n3 = floor((NgridZ-0)*(coor(2)-min(2))/(max(2)-min(2)));

            // Wrap if extends over the grid dimensions
            while(n1>=NgridX || n1<0)
                n1>=0 ? n1 %= NgridX : n1 = NgridX + n1%NgridX;
            while(n2>=NgridY || n2<0)
                n2>=0 ? n2 %= NgridY : n2 = NgridY + n2%NgridY;
            while(n3>=NgridZ || n3<0)
                n3>=0 ? n3 %= NgridZ : n3 = NgridZ + n3%NgridZ;

            // Assign to grid
            grid[n1][n2][n3].push_back(i);
        }
    }

    /*
    // Statistics:
    int min = 1e20,max = -1e20, cur, zero = 0, tot = 0;
    for(i=0;i<NgridX;++i)
        for(j=0;j<NgridY;++j)
            for(k=0;k<NgridZ;++k){
                cur = grid[i][j][k].size();
                if(cur<min) min = cur;
                if(cur>max) max = cur;
                if(cur==0) ++zero;
                tot += cur;
            }
    cout << "Number of cells: " << NgridX*NgridY*NgridZ << endl;
    cout << "Number of effective atoms: " << tot << endl;
    cout << "Minimal per cell: " << min << endl;
    cout << "Maximal per cell: " << max << endl;
    cout << "Number of empty cells: " << zero << endl << endl;
    */
}


// Search inside one selection
void Grid_searcher::do_search(Selection& sel, std::vector<Eigen::Vector2i>& bon,
                              std::vector<float>* dist_vec){
    int i,j,k,i1;
    int nlist_size;

    // Search
    bon.clear();
    if(dist_vec) dist_vec->clear();

    // Init visited cells array
    for(i=0;i<NgridX;++i)
        for(j=0;j<NgridY;++j)
            for(k=0;k<NgridZ;++k)
                visited[i][j][k] = false;

    // Searching
    for(i=0;i<NgridX;++i){
        for(j=0;j<NgridY;++j){
            for(k=0;k<NgridZ;++k){                
                // Search in central cell
                get_central_1(i,j,k, sel, bon, dist_vec);
                visited[i][j][k] = true;
                // Get neighbour list
                get_nlist(i,j,k);
                nlist_size = nlist.size();
                // Searh between this and neighbouring cells
                for(i1=0;i1<nlist_size;++i1){
                    //cout << nlist[i1].transpose() << endl;
                    if( !visited[nlist[i1](0)][nlist[i1](1)][nlist[i1](2)] )
                        get_side_1(i,j,k, nlist[i1](0),nlist[i1](1),nlist[i1](2),sel, bon, dist_vec);
                }
            }
        }
    }
    //cout << "Searcher has: " << bonds.size() << " pairs" << endl;
}

// Search between two selections
void Grid_searcher::do_search(Selection& sel1, Selection& sel2, std::vector<Eigen::Vector2i>& bon,
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


    for(i=0;i<NgridX;++i){
        for(j=0;j<NgridY;++j){
            for(k=0;k<NgridZ;++k){
                // Search in central cell
                get_central_2(i,j,k, sel1, sel2, bon, dist_vec);
                visited[i][j][k] = true;
                // Get neighbour list
                get_nlist(i,j,k);
                nlist_size = nlist.size();
                // Searh between this and neighbouring cells
                for(i1=0;i1<nlist_size;++i1){
                    //cout << nlist[i1].transpose()<<endl;
                    if( !visited[nlist[i1](0)][nlist[i1](1)][nlist[i1](2)] )
                        get_side_2(i,j,k, nlist[i1](0),nlist[i1](1),nlist[i1](2),
                                   sel1, sel2, bon, dist_vec);
                }

            }
        }
    }
    //cout << "Searcher has: " << bonds.size() << " pairs" << endl;
}

// Get periodic distance between two points for current box_dim
float Grid_searcher::periodic_distance(const Vector3f& p1, const Vector3f& p2){
    // For each dimension measure periodic distance
    Vector3f v = (p2-p1).array().abs();
    for(int i=0;i<3;++i)
        if(v(i)>0.5*box_dim(i)) v(i) = box_dim(i)-v(i);
    return v.norm();
}

// Search in central cell inside 1 selection
void Grid_searcher::get_central_1(int i1, int j1, int k1, Selection& sel,
                                    std::vector<Eigen::Vector2i>& bonds,
                                    std::vector<float>* dist_vec){
    int c1,c2;
    int n1 = grid1[i1][j1][k1].size();

    if(n1==0) return; //Nothing to do

    float d;
    Vector2i pair;
    Vector3f p1,p2;

    for(c1=0;c1<n1-1;++c1)
        for(c2=c1+1;c2<n1;++c2){
            if(!is_periodic)
                // Get non-periodic distance
                d = (sel.XYZ(grid1[i1][j1][k1][c1]) -
                     sel.XYZ(grid1[i1][j1][k1][c2])).norm();
            else
                // Get periodic distance
                d = periodic_distance(sel.XYZ(grid1[i1][j1][k1][c1]),
                                      sel.XYZ(grid1[i1][j1][k1][c2]));

            if(d<=cutoff){
                if(abs_index){
                    pair(0) = sel.Index(grid1[i1][j1][k1][c1]);
                    pair(1) = sel.Index(grid1[i1][j1][k1][c2]);
                } else {
                    pair(0) = grid1[i1][j1][k1][c1];
                    pair(1) = grid1[i1][j1][k1][c2];
                }

                bonds.push_back(pair); //Add bond
                if(dist_vec) dist_vec->push_back(d);
            }
        }
}

// Add bonds between two given cells of the grid for single selection
void Grid_searcher::get_side_1(int i1,int j1,int k1, int i2,int j2,int k2, Selection& sel,
                                std::vector<Eigen::Vector2i>& bonds,
                                std::vector<float>* dist_vec){

    int c1,c2;
    int n1 = grid1[i1][j1][k1].size();
    int n2 = grid1[i2][j2][k2].size();

    if(n1==0 || n2==0) return; //Nothing to do

    float d;
    Vector2i pair;

    for(c1=0;c1<n1;++c1)
        for(c2=0;c2<n2;++c2){
            if(!is_periodic)
                d = (sel.XYZ(grid1[i1][j1][k1][c1]) -
                     sel.XYZ(grid1[i2][j2][k2][c2])).norm();
            else
                d = periodic_distance(sel.XYZ(grid1[i1][j1][k1][c1]),
                                      sel.XYZ(grid1[i2][j2][k2][c2]));

            if(d<=cutoff){
                if(abs_index){
                    pair(0) = sel.Index(grid1[i1][j1][k1][c1]);
                    pair(1) = sel.Index(grid1[i2][j2][k2][c2]);
                } else {
                    pair(0) = grid1[i1][j1][k1][c1];
                    pair(1) = grid1[i2][j2][k2][c2];
                }
                bonds.push_back(pair);
                if(dist_vec) dist_vec->push_back(d);
            }
        }
}


// Get list of neighbouring cells. Depends on cutoff!
void Grid_searcher::get_nlist(int i,int j,int k){
    Vector3i coor;
    int c1,c2,c3;
    float x,y,z;

    nlist.clear();

    // Determine how many cells the cutoff span in directions X, Y and Z
    int spanX = ceil(cutoff/dX);
    int spanY = ceil(cutoff/dY);
    int spanZ = ceil(cutoff/dZ);    

    // In degenerate cases of one cell for some dimension we need to
    // suppress this dimension by leaving only central cell
    if(NgridX==1) spanX = 0;
    if(NgridY==1) spanY = 0;
    if(NgridZ==1) spanZ = 0;

    //cout << "Spans " << spanX << " " << spanY << " " << spanZ << endl;

    if(is_periodic)
        // Periodic variant
        for(c1=-spanX;c1<=spanX;++c1)
            for(c2=-spanY;c2<=spanY;++c2)
                for(c3=-spanZ;c3<=spanZ;++c3){
                    coor(0) = i+c1;
                    coor(1) = j+c2;
                    coor(2) = k+c3;                                       

                    //Exclude central cell
                    if(coor(0) == i && coor(1) == j && coor(2) == k ) continue;

                    //cout << "row " << coor.transpose() << endl;

                    // Periodic variant: wrap cells around
                    while(coor(0)>=NgridX || coor(0)<0)
                        coor(0)>=0 ? coor(0) %= NgridX : coor(0) = NgridX + coor(0)%NgridX;
                    while(coor(1)>=NgridY || coor(1)<0)
                        coor(1)>=0 ? coor(1) %= NgridY : coor(1) = NgridY + coor(1)%NgridY;
                    while(coor(2)>=NgridZ || coor(2)<0)
                        coor(2)>=0 ? coor(2) %= NgridZ : coor(2) = NgridZ + coor(2)%NgridZ;

                    //cout << "wrapped " << coor.transpose() << endl;
                    // This may be a corner cell, which is not included in fact,
                    // but currently we ignore this
                    nlist.push_back(coor); // Add cell
                }
    else
        // Non-periodic variant
        for(c1=-spanX;c1<=spanX;++c1)
            for(c2=-spanY;c2<=spanY;++c2)
                for(c3=-spanZ;c3<=spanZ;++c3){
                    coor(0) = i+c1;
                    coor(1) = j+c2;
                    coor(2) = k+c3;

                    //Bounds check
                    if(coor(0)<0 || coor(1)<0 || coor(2)<0 ||
                        coor(0)>=NgridX || coor(1)>=NgridY || coor(2)>=NgridZ) continue;

                    if(coor(0) == i && coor(1) == j && coor(2) == k ) continue; //Exclude central cell
                    // This may be a corner cell, which is not included in fact,
                    // but currently we ignore this
                    nlist.push_back(coor); // Add cell
                }
}


// Add bonds inside one cell for two selections
void Grid_searcher::get_central_2(int i1,int j1,int k1, Selection& sel1, Selection& sel2,
                                    std::vector<Eigen::Vector2i>& bonds,
                                    std::vector<float>* dist_vec){
    int c1,c2;
    int n1 = grid1[i1][j1][k1].size();
    int n2 = grid2[i1][j1][k1].size();
    if(n1==0 || n2==0) return; //Nothing to do

    float d;
    VectorXi pair(2);

    for(c1=0;c1<n1;++c1)
        for(c2=0;c2<n2;++c2){
            if(!is_periodic)
                d = (sel1.XYZ(grid1[i1][j1][k1][c1]) - sel2.XYZ(grid2[i1][j1][k1][c2])).norm();
            else
                d = periodic_distance(sel1.XYZ(grid1[i1][j1][k1][c1]),
                                      sel2.XYZ(grid2[i1][j1][k1][c2]));

            if(d<=cutoff){
                if(abs_index){
                    pair(0) = sel1.Index(grid1[i1][j1][k1][c1]);
                    pair(1) = sel2.Index(grid2[i1][j1][k1][c2]);
                } else {
                    pair(0) = grid1[i1][j1][k1][c1];
                    pair(1) = grid2[i1][j1][k1][c2];
                }
                bonds.push_back(pair);
                if(dist_vec) dist_vec->push_back(d);
            }
        }
}

// Add contacts between two given cells of the grid for two selections
void Grid_searcher::get_side_2(int i1,int j1,int k1, int i2,int j2,int k2,
                                Selection& sel1, Selection& sel2,
                                std::vector<Eigen::Vector2i>& bonds,
                                std::vector<float>* dist_vec){
    int c1,c2;
    float d;
    Vector2i pair;

    // First phase. Search for contacts between sel1 in cell1 and sel2 in cell2
    int n1 = grid1[i1][j1][k1].size();
    int n2 = grid2[i2][j2][k2].size();

    if(n1*n2!=0){ // If both cells are not empty
        for(c1=0;c1<n1;++c1)
            for(c2=0;c2<n2;++c2){

                if(!is_periodic)
                    d = (sel1.XYZ(grid1[i1][j1][k1][c1]) - sel2.XYZ(grid2[i2][j2][k2][c2])).norm();
                else {
                    d = periodic_distance(sel1.XYZ(grid1[i1][j1][k1][c1]),
                                          sel2.XYZ(grid2[i2][j2][k2][c2]));
                }

                if(d<=cutoff){
                    if(abs_index){
                        pair(0) = sel1.Index(grid1[i1][j1][k1][c1]);
                        pair(1) = sel2.Index(grid2[i2][j2][k2][c2]);
                    } else {
                        pair(0) = grid1[i1][j1][k1][c1];
                        pair(1) = grid2[i2][j2][k2][c2];
                    }
                    bonds.push_back(pair);
                    if(dist_vec) dist_vec->push_back(d);
                }
            }
    }

    // Second phase. Search for contacts between sel2 in cell1 and sel1 in cell2
    n1 = grid2[i1][j1][k1].size();
    n2 = grid1[i2][j2][k2].size();

    if(n1*n2!=0){ // If cells are not empty
        for(c1=0;c1<n1;++c1)
            for(c2=0;c2<n2;++c2){
                if(!is_periodic)
                    d = (sel2.XYZ(grid2[i1][j1][k1][c1]) - sel1.XYZ(grid1[i2][j2][k2][c2])).norm();
                else {
                    d = periodic_distance(sel2.XYZ(grid2[i1][j1][k1][c1]),
                                          sel1.XYZ(grid1[i2][j2][k2][c2]));
                }

                if(d<=cutoff){
                    if(abs_index){
                        pair(1) = sel2.Index(grid2[i1][j1][k1][c1]); //ordered pair!
                        pair(0) = sel1.Index(grid1[i2][j2][k2][c2]);
                    } else {
                        pair(1) = grid2[i1][j1][k1][c1]; //ordered pair!
                        pair(0) = grid1[i2][j2][k2][c2];
                    }
                    bonds.push_back(pair);
                    if(dist_vec) dist_vec->push_back(d);
                }
            }
    }
}
