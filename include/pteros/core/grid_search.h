/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
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

    /** Implements grid search algorithm for finding atoms, located within
      certain distance from their centers.
      Returns pairs of "contacting" atoms
       as a plain list of index pairs
    */
    class Grid_searcher {
        public:
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

            /// Default constructor
            Grid_searcher();

            /// Assign atoms from given selection to grid
            /// Pointer to this selection is stored internally
            void assign_to_grid(float d, Selection& sel,
                                bool absolute_index = false,
                                bool periodic = false);

            /// Search contacts within distance from given point in space
            /// Existing grid set in assign_to_grid() is used
            /// Returns the list of atoms in contact with given point
            void search_within(Eigen::Vector3f& coor,
                              std::vector<int>& bon);

            /// Search within some distance from target selection
            /// sel must be the same as used in assign_to_grid() to give meaningful results
            void search_within(Selection& target,
                               std::vector<int>& bon,
                               bool include_self=true);

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

            /*
            // Coordinate accessors
            // Could be overloaded to provede access to coordinates
            // stored not in the system itself but somewhere
            // Passed selection is used to extract absolute index in any case.
            inline virtual float getX(const int i, Selection& sel){
                return sel.X(i);
            }

            inline virtual float getY(const int i, Selection& sel){
                return sel.Y(i);
            }

            inline virtual float getZ(const int i, Selection& sel){
                return sel.Z(i);
            }

            inline virtual Eigen::Vector3f getXYZ(const int i, Selection& sel){
                return sel.XYZ(i);
            }
            */
    };

}

#endif // GRID_SEARCH_H_INCLUDED
