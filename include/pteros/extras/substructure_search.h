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

#include <vector>
#include "pteros/core/selection.h"

namespace pteros {

/// Finds molecular symmetry using OpenBabel based on bonding pattern.
/// Returns groups of atoms which are topologically equivalent (local indexes are returned)
std::vector<std::vector<int>> find_equivalent_atoms(const Selection& sel, int x_memory);

/// Find mapping between source and query selection.
/// Size of query have to be <= size of source.
/// Returned array shows mapping from query to source such as
/// atom i in query corresponds to result[i] in source
/// Indexes are local to both selections.
/// @param find_all - if true returns all unique mappings. If false returns only the first mapping.
std::vector<std::vector<int>> find_substructures(const Selection& source, const Selection& query, bool find_all=false);

/// Make target topologically equivalent to templ and return it as a new system.
/// Resulting system has coordinates of target but the sequence of atoms of template.
/// Missed hydrogens are added to target prior to rearrangement if needed.
/// Typical usage is to make docking poses compatible with all-atom topology for MD.
System make_equivalent_to_template(const Selection& target, const Selection& templ);

} // namespace





