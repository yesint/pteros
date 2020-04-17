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



#ifndef DISTANCE_SEARCH_WITHIN_BASE_H_INCLUDED
#define DISTANCE_SEARCH_WITHIN_BASE_H_INCLUDED

#include "distance_search_base.h"
#include "atomic_wrapper.h"

namespace pteros {       

class Distance_search_within_base: public Distance_search_base {
protected:
    // Array of atomic bools for used source points
    std::vector<atomic_wrapper<bool>> used;
    void do_search(int sel_size);
    void do_part(int dim, int _b, int _e);
    // Pointer to source selection
    Selection* src_ptr;
    void search_in_pair_of_cells(int sx, int sy, int sz, int tx, int ty, int tz, bool is_periodic);
    void used_to_result(std::vector<int>& res, bool include_self,
                        const Selection &src, const Selection &target);
};

}

#endif


