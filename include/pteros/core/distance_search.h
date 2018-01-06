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


/// @file

#ifndef DISTANCE_SEARCH_INCLUDED
#define DISTANCE_SEARCH_INCLUDED

#include "pteros/core/selection.h"
#include "pteros/core/distance_search_within.h"

namespace pteros {       

/// Search contacts within single selection
/// Optionally returns distances for each pair.
void search_contacts(float d,
                     const Selection& sel,
                     std::vector<Eigen::Vector2i> &pairs,
                     bool absolute_index = false,
                     bool periodic = false,
                     std::vector<float> *dist_vec = nullptr);

/// Search contacts between two selections
void search_contacts(float d,
                     const Selection& sel1,
                     const Selection& sel2,
                     std::vector<Eigen::Vector2i>& pairs,
                     bool absolute_index = false,
                     bool periodic = false,
                     std::vector<float>* dist_vec = nullptr);

/// Search atoms from source selection around the traget selection
/// Returns absolute indexes only!
void search_within(float d,
                   const Selection& src,
                   const Selection& target,
                   std::vector<int> &res,
                   bool include_self=true,
                   bool periodic = false);

}

#endif

