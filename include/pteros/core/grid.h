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


#pragma once

#include <Eigen/Core>
#include <vector>
#include <deque>
#include "pteros/core/selection.h"
#include <unsupported/Eigen/CXX11/Tensor>

namespace pteros {    

    class GridCell {
    public:
        void add_point(int ind, Vector3f_const_ref crd);
        void clear();
        int get_index(int i) const {return indexes[i];}
        Eigen::Vector3f get_coord(int i) const {return coords[i];}
        size_t size() const {return indexes.size();}

    private:
        std::vector<int> indexes;
        std::vector<Eigen::Vector3f> coords;
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
        GridCell& cell(int i, int j, int k){ return data(i,j,k); }
        GridCell& cell(Vector3i_const_ref ind){ return data(ind(0),ind(1),ind(2)); }
        const GridCell& cell(Vector3i_const_ref ind) const { return data(ind(0),ind(1),ind(2)); }

        /// Non-periodic populate
        void populate(const Selection& sel,bool abs_index = false);

        void populate(const Selection& sel,
                      Vector3f_const_ref min,
                      Vector3f_const_ref max,
                      bool abs_index);

        /// Periodic populate
        void populate_periodic(const Selection& sel,
                               Vector3i_const_ref pbc_dims = fullPBC,
                               bool abs_index = false);

        void populate_periodic(const Selection& sel,
                               const PeriodicBox& box,
                               Vector3i_const_ref pbc_dims = fullPBC,
                               bool abs_index = false);
    private:
        Eigen::Tensor<GridCell,3> data;
    };

}




