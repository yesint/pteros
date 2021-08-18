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

#include <string>
#include <vector>

namespace pteros {

/// Macro definitions for selections. Each macro is expanded during
/// evaluation of selection.
static const std::vector<std::string> selection_macro {
    "protein", "(resname ABU ACE AIB ALA ARG ARGN ASN ASN1 ASP ASP1 ASPH CYS CYS1 CYS2 CYSH DALA GLN GLU GLUH GLY HIS HIS1 HISA HISB HISH HSD HSE HSP HYP ILE LEU LYS LYSH MELEU MET MEVAL NAC NH2 PHE PHEH PHEU PHL PRO SER THR TRP TRPH TRPU TYR TYRH TYRU VAL PGLU)",
    "backbone", "(name N CA C and resname ABU ACE AIB ALA ARG ARGN ASN ASN1 ASP ASP1 ASPH CYS CYS1 CYS2 CYSH DALA GLN GLU GLUH GLY HIS HIS1 HISA HISB HISH HSD HSE HSP HYP ILE LEU LYS LYSH MELEU MET MEVAL NAC NH2 PHE PHEH PHEU PHL PRO SER THR TRP TRPH TRPU TYR TYRH TYRU VAL PGLU)",
    "acidic", "(resname ASP GLU)",
    "cyclic", "(resname HIS PHE PRO TRP TYR)",
    "aromatic", "(resname HIS PHE TRP TYR)",
    "basic", "(resname ARG HIS LYS HSP)",
    "buried", "(resname ALA LEU VAL ILE PHE CYS MET TRP)",
    "charged", "(resname ARG HIS LYS HSP ASP GLU)",
    "hydrophobic", "(resname ALA LEU VAL ILE PRO PHE MET TRP)",
    "water", "(resname HOH SOL TIP3)",
    "hydrogen", "(name 'H.*')",
    "heavy", "(not name 'H.*')",
    "noh", "(not name 'H.*')",
    "DNA", "(resname A T G C U)",
    "now", "(not resname HOH SOL TIP3)",
};

} // end of namespace pteros




