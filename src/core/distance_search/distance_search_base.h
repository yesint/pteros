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

#pragma once

#include <Eigen/Core>
#include <vector>
#include "pteros/core/periodic_box.h"
#include "pteros/core/grid.h"

namespace pteros {

    class Distance_search_base {    
    protected:
        // Min and max of the bounding box (for non-periodic case)
        Eigen::Vector3f min,max;
        // Current periodic box (for periodic case)
        Periodic_box box;
        // Grid dimensions
        int NgridX, NgridY, NgridZ;
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

        // Search plan
        void make_search_plan(std::vector<Eigen::Matrix<int,3,2>>& plan);
        // Periodic grid size
        void set_grid_size(const Eigen::Vector3f& min,
                           const Eigen::Vector3f& max);
        // Non-periodic grid size
        void set_grid_size(const Periodic_box& box);
        // Create single grid
        void create_grid(const Selection &sel);
        // Create two grids
        void create_grids(const Selection &sel1, const Selection &sel2);

    };

}


