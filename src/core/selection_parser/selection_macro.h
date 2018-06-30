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

#ifndef SELECTION_MACRO_H
#define SELECTION_MACRO_H

#include <string>
#include <vector>

namespace pteros {

/// Macro definitions for selections. Each macro is expanded during
/// evaluation of selection.
static const std::vector<std::string> macro {
    "protein", "(resname ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR LYSH HISA HISB CYSH)",
    "backbone", "(name C CA O N)",
    "acidic", "(resname ASP GLU)",
    "cyclic", "(resname HIS PHE PRO TRP TYR)",
    "aromatic", "(resname HIS PHE TRP TYR)",
    "basic", "(resname ARG HIS LYS HSP)",
    "buried", "(resname ALA LEU VAL ILE PHE CYS MET TRP)",
    "charged", "(resname ARG HIS LYS HSP ASP GLU)",
    "hydrophobic", "(resname ALA LEU VAL ILE PRO PHE MET TRP)",
    "water", "(resname HOH SOL)"
};

} // end of namespace pteros

#endif

