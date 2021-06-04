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

#include "pteros/core/system.h"

namespace pteros {

// Mol_file is a friend of System and can access it's internals
// but derived *_file classes are not friends.
// In order to access internals of the System we define special access class
class SystemBuilder {
public:
    SystemBuilder(System& s): sys(&s) {}
    SystemBuilder(System* s): sys(s) {}
    // When destroyed builer calls assign_resindex() and duing other preparations
    ~SystemBuilder();

    void allocate_atoms(int n);
    void set_atom(int i, const Atom& at);
    Atom& atom(int i);
    void add_atom(const Atom& at);
private:
    System* sys;
};

} // namespace





