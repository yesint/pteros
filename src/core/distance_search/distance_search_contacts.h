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


#ifndef DISTANCE_SEARCH_CONTACTS_H_INCLUDED
#define DISTANCE_SEARCH_CONTACTS_H_INCLUDED

#include "distance_search_base.h"

namespace pteros {       

class Distance_search_contacts: public Distance_search_base {
public:
protected:
    // wisited array
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

