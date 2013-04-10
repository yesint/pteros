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
#ifndef SELECTION_MACRO_H
#define SELECTION_MACRO_H

namespace pteros {

/// Macro definitions for selections. Each macro is expanded during
/// evaluation of selection.
static const char* macro[] = {
    "backbone", "(name C CA O N)",
    "acidic", "(resname ASP GLU)",
    "cyclic", "(resname HIS PHE PRO TRP TYR)",
    "aromatic", "(resname HIS PHE TRP TYR)",
    "basic", "(resname ARG HIS LYS HSP)",
    "buried", "(resname ALA LEU VAL ILE PHE CYS MET TRP)",
    "charged", "(resname ARG HIS LYS HSP ASP GLU)",
    "hydrophobic", "(resname ALA LEU VAL ILE PRO PHE MET TRP)"
};
/// Number of macro definitions
static const int Nmacro = 8;

} // end of namespace pteros

#endif
