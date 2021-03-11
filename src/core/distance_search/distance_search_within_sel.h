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

#include "distance_search_within_base.h"

namespace pteros {       


class DistanceSearchWithinSel: public DistanceSearchWithinBase {
public:
    /// Constuctor for very fast immediate search of atoms from src,
    /// which are within given distance from the atoms of target.
    /// Used in internal parsing of within selections.
    /// \warning Returns absolute indexes only!
    /// \warning Result is not sorted!
    DistanceSearchWithinSel(float d,
                           const Selection& src,
                           const Selection& target,
                           std::vector<int> &res,
                           bool include_self=true,
                           Vector3i_const_ref pbc = noPBC);
};

}




