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

#include "pteros/core/selection.h"
#include "boost/multi_array.hpp"

#include <Eigen/Core>
#include <vector>

namespace pteros {

    typedef boost::multi_array<std::vector<int>,3> Grid_t;

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
    \endcode
    - Sorting the atoms from given selection into the cells of
    3D grid with given dimensions. Useful for producing volumetric datasets or
    various 3D histrograms (for examples density or the residence time maps).
    Typical usage:
    \code
    Grid_searcher g;
    g.create_custom_grid(NX,NY,NZ);
    g.fill_custom_grid(sel,true);
    for(i=0;i<NX;++i)
        for(j=0;j<NY;++j)
            for(k=0;k<NZ;++k)
                cout << g.cell_of_custom_grid(i,j,k) << endl;
    \endcode
    */

    class Grid_searcher {
        public:
            /// Default constructor
            Grid_searcher();

            /// @name Search for atoms located within the cutoff distance from each other
            /// @{

            /// Constructor for instant searching inside a selection
            Grid_searcher(float d, Selection& sel,
                            std::vector<Eigen::Vector2i>& bon,                            
                            bool absolute_index = false,
                            bool periodic = false,
                            std::vector<float>* dist_vec = NULL);

            /// Constructor for instant searching between two selections
            Grid_searcher(float d, Selection& sel1, Selection& sel2,
                            std::vector<Eigen::Vector2i>& bon,
                            bool absolute_index = false,
                            bool periodic = false,
                            std::vector<float>* dist_vec = NULL);
            /// @}

            /// @name Search for atoms within the given distance from point or selection
            /// @{

            /// Assign atoms from given selection to grid
            /// Pointer to this selection is stored internally
            void assign_to_grid(float d, Selection& sel,
                                bool absolute_index = false,
                                bool periodic = false);

            /// Search atoms within given distance from point in space
            /// Existing grid set in assign_to_grid() is used
            /// Returns the list of atoms in the vicinity of given point
            void search_within(Eigen::Vector3f& coor,
                              std::vector<int>& bon);

            /// Search atoms within given distance from target selection
            /// target must be the subset of selectionm which was used in assign_to_grid()
            /// to give meaningful results
            void search_within(Selection& target,
                               std::vector<int>& bon,
                               bool include_self=true);
            /// @}

            /// @name Assigning the atoms from given selection to the custom periodic grid
            /// @{

            /// Creates custom periodic grid with given dimensions
            void create_custom_grid(int nX, int nY, int nZ);
            /// Populate custom grid created by create_custom_grid() from given selection
            void fill_custom_grid(Selection sel,
                                  bool absolute_index = false);
            /// Read/write acces to the cells of custom grid
            std::vector<int>& cell_of_custom_grid(int x, int y, int z);

            /// @}


        protected:
            // Create one grid from single selection
            void create_grid(Grid_t& grid, Selection& sel);
            // Create two grids from two selections
            void create_grid2(Selection& sel1, Selection& sel2);

            void populate_grid(Grid_t& grid, Selection& sel);

            /// Search function for contacts inside one group
            void do_search(Selection& sel, std::vector<Eigen::Vector2i>& bon,
                           std::vector<float>* dist_vec);
            /// Search function for contacts between two groups
            void do_search(Selection& sel1, Selection& sel2, std::vector<Eigen::Vector2i>& bon,
                           std::vector<float>* dist_vec);
            // Min and max of the bounding box
            Eigen::Vector3f min,max;
            // Grid dimensiond
            int NgridX, NgridY, NgridZ;
            // Grid periodic steps
            float dX, dY, dZ;
            // Grids
            Grid_t grid1, grid2;
            boost::multi_array<bool, 3> visited;
            // Neighbour list for cells
            std::vector<Eigen::Vector3i> nlist;
            // Cut-off
            float cutoff;

            // Pointer to selection used in within searching
            Selection* p_sel;

            // Basis conversion stuff for triclinic boxes
            Eigen::Matrix3f inv_basis_matr;
            // Current periodic box sizes
            Eigen::Vector3f box_dim;
            void make_inv_matr(const Eigen::Matrix3f& box);

            // If true absolute index rather then selection index is returned in the bond list
            bool abs_index;
            bool is_periodic;
            bool is_triclinic;

            // Periodic distance between two points
            inline float periodic_distance(const Eigen::Vector3f& p1, const Eigen::Vector3f& p2);

            void set_grid_size(const Eigen::Vector3f& min, const Eigen::Vector3f& max, int Natoms);

            void get_nlist(int i,int j,int k);            

            void get_central_1(int i1, int j1, int k1, Selection& sel,
                                std::vector<Eigen::Vector2i>& bonds,
                                std::vector<float>* dist_vec);
            void get_side_1(int i1,int j1,int k1, int i2,int j2,int k2, Selection& sel,
                                std::vector<Eigen::Vector2i>& bonds,
                                std::vector<float>* dist_vec);
            void get_central_2(int i1, int j1, int k1, Selection& sel1, Selection& sel2,
                                std::vector<Eigen::Vector2i>& bonds,
                                std::vector<float>* dist_vec);
            void get_side_2(int i1,int j1,int k1, int i2,int j2,int k2,
                                Selection& sel1,
                                Selection& sel2,
                                std::vector<Eigen::Vector2i>& bonds,
                                std::vector<float>* dist_vec);
    };

}

#endif // GRID_SEARCH_H_INCLUDED
