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

#ifndef GRID_SEARCH_H_INCLUDED
#define GRID_SEARCH_H_INCLUDED

#include <Eigen/Core>
#include <vector>
#include "pteros/core/selection.h"
#include "pteros/core/atomic_wrapper.h"
#include "pteros/core/grid.h"

namespace pteros {       

    struct Nlist_t {
        std::vector<Eigen::Vector3i> data;
        std::vector<bool> wrapped;

        void clear();
        void append(Vector3i_const_ref coor, bool wrap = false);
    };

    /** @brief Implements grid search algorithm
    Grid_searcher class subdivides the volume of the system into number of
    cells and computes which atoms appear in each cell. After that the pairs of
    atoms, which reside within the given cut-off distance from each other
    could be computed rapidly. Grid search method is much faster than simple brute
    force search.
    Grid_searcher class could be used in three different scenarios:

    - Search for atoms, which are located within the cutoff
    inside given selection or between two given selections. Typical usage:
    \code
    vector<Vector2i> bonds;
    Grid_searcher(0.2, sel, bonds, true);
    Grid_searcher(0.2, sel1, sel2, bonds, true);
    \endcode
    - Search for atoms, which are located within the given distance from some
    arbitrary point in space or from any atom of the given selection.
    This scenario is also used internally when processing "within" keywords of selection syntax.
    Typical usage:
    \code
    vector<int> atoms;
    Selection all(system,"all");
    Selection res_of_interest(system,"resid 100");

    Grid_searcher g;
    g.assign_to_grid(0.2, all, true, true);
    g.search_within(res_of_interest, atoms);

    // The same using constructor with immediate searching:
    Grid_searcher(0.2, all, res_of_interest, atoms, false, true, true);
    \endcode    
    */

    class Grid_searcher {
        public:            
            /// Default constructor
            Grid_searcher();

            /// @name Search for atoms located within the cutoff distance from each other
            /// @{

            /// Constructor for instant searching inside a selection
            Grid_searcher(float d, const Selection& sel,
                            std::vector<Eigen::Vector2i>& bon,                            
                            bool absolute_index = false,
                            bool periodic = false,
                            std::vector<float>* dist_vec = NULL);

            /// Constructor for instant searching between two selections
            Grid_searcher(float d, const Selection& sel1, const Selection& sel2,
                            std::vector<Eigen::Vector2i>& bon,
                            bool absolute_index = false,
                            bool periodic = false,
                            std::vector<float>* dist_vec = NULL);
            /// @}

            /// @name Search for atoms within the given distance from point or selection
            /// @{

            /// Assign atoms from given selection to grid
            /// Pointer to this selection is stored internally
            void assign_to_grid(float d, const Selection& sel,
                                bool absolute_index = false,
                                bool periodic = false);

            /// Search atoms within given distance from point in space
            /// Existing grid set in assign_to_grid() is used
            /// Returns the list of atoms in the vicinity of given point
            void search_within(Vector3f_const_ref coord,
                              std::vector<int>& bon);

            /// Search atoms within given distance from target selection
            /// target must be the subset of selection, which was used in assign_to_grid()
            /// to give meaningful results
            void search_within(const Selection& target,
                               std::vector<int>& bon,
                               bool include_self=true);

            /// Constuctor for very fast immediate search of atoms from src,
            /// which are within given distance from the atoms of target.
            /// Used in internal parsing of within selections.
            /// \warning Returns absolute indexes only!
            Grid_searcher(  float d,
                            const Selection& src,
                            const Selection& target,
                            std::vector<int>& bon,
                            bool include_self=true,                            
                            bool periodic = false);
            /// @}


        protected:

            // Create one grid from single selection            
            void create_grid(Grid& grid, const Selection& sel);
            // Create two grids from two selections            
            void create_grid2(const Selection& sel1, const Selection& sel2);            

            /// Search function for contacts inside one group
            void do_search1(std::vector<Eigen::Vector2i>& bon,
                           std::vector<float>* dist_vec);
            /// Search function for contacts between two groups
            void do_search2(std::vector<Eigen::Vector2i>& bon,
                           std::vector<float>* dist_vec);

            void do_search_within(std::vector<int>& bon, const Selection& src);

            void do_part1(int dim, int _b, int _e,
                          std::vector<Eigen::Vector2i>& bon,
                          std::vector<float>* dist_vec);

            void do_part2(int dim, int _b, int _e,
                          std::vector<Eigen::Vector2i>& bon,
                          std::vector<float>* dist_vec);

            void do_part_within(int dim, int _b, int _e,
                                std::vector<atomic_wrapper<bool>>& used);


            // Min and max of the bounding box
            Eigen::Vector3f min,max;
            // Grid dimensions
            int NgridX, NgridY, NgridZ;

            // Grids with coordinate pointers
            Grid grid1,grid2;

            boost::multi_array<bool, 3> visited;

            // Cut-off
            float cutoff;

            // Pointer to selection used in within searching
            Selection* p_sel;

            // Current periodic box
            Periodic_box box;

            // If true absolute index rather then selection index is returned in the bond list
            bool abs_index;
            bool is_periodic;            

            void set_grid_size(const Eigen::Vector3f& min, const Eigen::Vector3f& max,
                               int Natoms, const Periodic_box& box);

            void get_nlist(int i, int j, int k, Nlist_t &nlist);
    };

}

#endif // GRID_SEARCH_H_INCLUDED
