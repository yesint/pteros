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

#ifndef DISTANCE_SEARCH_BASE_H_INCLUDED
#define DISTANCE_SEARCH_BASE_H_INCLUDED

#include <Eigen/Core>
#include <vector>
#include "pteros/core/periodic_box.h"
#include "pteros/core/grid.h"

namespace pteros {       

    struct Nlist_t {
        std::vector<Eigen::Vector3i> data;
        std::vector<bool> wrapped;

        void clear();
        void append(Vector3i_const_ref coor, bool wrap = false);
    };


    class Distance_search_base {
    public:
    protected:
        // Min and max of the bounding box (for non-periodic case)
        Eigen::Vector3f min,max;
        // Current periodic box (for periodic case)
        Periodic_box box;
        // Grid dimensions
        int NgridX, NgridY, NgridZ;
        // Grids with coordinate pointers
        Grid grid1,grid2;
        // Cut-off
        float cutoff;
        // If true absolute index rather then selection index is returned in the bond list
        bool abs_index;
        // Periodicity
        bool is_periodic;

        void set_grid_size(const Eigen::Vector3f& min, const Eigen::Vector3f& max,
                           int Natoms, const Periodic_box& box);


        void get_nlist(int i, int j, int k, Nlist_t &nlist);
    };

}

#endif // GRID_SEARCH_H_INCLUDED
