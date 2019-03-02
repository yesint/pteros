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

#pragma once

#include <Eigen/Core>
#include <vector>
#include <deque>
#include "pteros/core/selection.h"

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

namespace pteros {    

    struct Grid_element {
        int index;
        Eigen::Vector3f* coor_ptr;
        Grid_element(int i, Eigen::Vector3f* ptr): index(i),coor_ptr(ptr) {}
    };

    /**
    Sorting the atoms from given selection into the cells of
    3D grid with given dimensions. Useful for producing volumetric datasets or
    various 3D histrograms (for examples density or the residence time maps).
    \code
    // Create 100x100x100 grid
    Grid g(100,100,100);
    // Populate it from given selection
    // in periodic manner
    g.populate_periodic(sel);

    // Print number of atoms in the grid cells
    for(int i=0;i< 100;++i){
        for(int j=0;j< 100;++j)
            for(int k=0;k< 100;++k)
                cout << g.cell(i,j,k).size() << endl;
    \endcode
     */
    class Grid {
    public:
        Grid(){}
        Grid(int X, int Y, int Z){ resize(X,Y,Z); }
        virtual ~Grid(){}

        void clear();
        void resize(int X, int Y, int Z);
        std::vector<Grid_element>& cell(int i, int j, int k){ return data[i][j][k]; }

        /// Non-periodic populate
        void populate(const Selection& sel,bool abs_index = false);

        void populate(const Selection& sel,
                      Vector3f_const_ref min,
                      Vector3f_const_ref max,
                      bool abs_index);

        /// Periodic populate
        void populate_periodic(const Selection& sel,bool abs_index = false);

        void populate_periodic(const Selection& sel,
                      const Periodic_box& box,
                      bool abs_index);
    private:        
        boost::multi_array<std::vector<Grid_element>,3> data;
        // Array of atomic coordinates, which have to be wrapped if periodic.
        // This is in order not to touch real coordinates of atoms and improve speed.
        std::deque<Eigen::Vector3f> wrapped_atoms;        
    };

}

