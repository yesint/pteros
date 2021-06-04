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
#include "pteros/core/periodic_box.h"
#include "pteros/core/grid.h"

namespace pteros {

    struct PlannedPair {
        Eigen::Vector3i c1;
        Eigen::Vector3i c2;
        Eigen::Vector3i wrapped;
    };


    class DistanceSearchBase {
    protected:
        // Min and max of the bounding box (for non-periodic case)
        Eigen::Vector3f min,max;
        // Current periodic box (for periodic case)
        PeriodicBox box;
        // Grid dimensions
        Eigen::Vector3i Ngrid;
        // Grids with coordinates
        Grid grid1,grid2;
        // Cut-off
        float cutoff;
        // If true absolute index rather then selection index is returned in the bond list
        bool abs_index;        
        // Periodic dimensions
        Eigen::Vector3i periodic_dims;
        // Is periodicity required?
        bool is_periodic;

        // Periodic grid size
        void set_grid_size(const Eigen::Vector3f& min,
                           const Eigen::Vector3f& max);
        // Non-periodic grid size
        void set_grid_size(const PeriodicBox& box);
        // Create single grid
        void create_grid(const Selection &sel);
        // Create two grids
        void create_grids(const Selection &sel1, const Selection &sel2);

        // Neighbours stencil
        const std::vector<Eigen::Vector3i> stencil = {
            // Center
            {0,0,0},{0,0,0},
            // Edges
            {0,0,0},{1,0,0}, //X
            {0,0,0},{0,1,0}, //Y
            {0,0,0},{0,0,1}, //Z
            // Face angles
            {0,0,0},{1,1,0}, //XY
            {0,0,0},{1,0,1}, //XZ
            {0,0,0},{0,1,1}, //YZ
            // Far angle
            {0,0,0},{1,1,1}, //XYZ
            // Face-diagonals
            {1,0,0},{0,1,0}, // XY
            {1,0,0},{0,0,1}, // XZ
            {0,1,0},{0,0,1}, // YZ
            // Cross-diagonals
            {1,1,0},{0,0,1}, // XY-Z
            {1,0,1},{0,1,0}, // XZ-Y
            {0,1,1},{1,0,0}, // YZ-X
        };

        Eigen::Vector3i index_to_pos(int i);
        bool process_neighbour_pair(PlannedPair &pair);
    };

}




