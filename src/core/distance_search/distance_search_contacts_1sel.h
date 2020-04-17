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



#ifndef DISTANCE_SEARCH_CONTACTS_1SEL_H_INCLUDED
#define DISTANCE_SEARCH_CONTACTS_1SEL_H_INCLUDED

#include "distance_search_contacts.h"

namespace pteros {       

class Distance_search_contacts_1sel: public Distance_search_contacts {
public:
    Distance_search_contacts_1sel(float d, const Selection& sel,
                                  std::vector<Eigen::Vector2i> &pairs,
                                  bool absolute_index = false,
                                  bool periodic = false,
                                  std::vector<float> *dist_vec = nullptr);
protected:
    void create_grid(const Selection &sel);

    void do_part(int dim, int _b, int _e,
                 std::deque<Eigen::Vector2i>& bon,
                 std::deque<float>* dist_vec) override;

    void search_in_cell(int x, int y, int z,
                        std::deque<Eigen::Vector2i> &bon,
                        std::deque<float> *dist_vec,
                        bool is_periodic);
};

}

#endif


