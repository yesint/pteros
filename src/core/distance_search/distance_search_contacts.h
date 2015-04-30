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

#ifndef DISTANCE_SEARCH_CONTACTS_H_INCLUDED
#define DISTANCE_SEARCH_CONTACTS_H_INCLUDED

#include "distance_search_base.h"

namespace pteros {       

class Distance_search_contacts: public Distance_search_base {
public:
protected:
    // withited array
    boost::multi_array<bool, 3> visited;
    // Pointers for final results
    std::vector<Eigen::Vector2i>* result_pairs;
    std::vector<float>* result_distances;

    void do_search();
    virtual void do_part(int dim, int _b, int _e,
                         std::deque<Eigen::Vector2i>& bon,
                         std::deque<float>* dist_vec) = 0;

    void search_in_pair_of_cells(int x1, int y1, int z1,
                                 int x2, int y2, int z2,
                                 Grid &grid1,
                                 Grid &grid2,
                                 std::deque<Eigen::Vector2i> &bon,
                                 std::deque<float> *dist_vec,
                                 bool is_periodic);
};

}

#endif
