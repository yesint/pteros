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

#include "distance_search_base.h"

namespace pteros {       

// In this class searching returns only the list of atoms from the source seletion.
// No distances are computed.
class DistanceSearchWithinBase: public DistanceSearchBase {
protected:
    // Implements logic for calling search_between_cells() or search_inside_cell()
    // with correct grids in derived classes

    // Pointer to search results
    std::unordered_set<int> result;

    void do_search();

    void compute_chunk(int b, int e, std::unordered_set<int> &res_buf);

private:
    void search_between_cells(Vector3i_const_ref c1,
                              Vector3i_const_ref c2,
                              Vector3i_const_ref wrapped,
                              std::unordered_set<int> &res_buffer);

    void search_planned_pair(const PlannedPair& pair,
                             std::unordered_set<int> &res_buffer);
};

}




