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


#include "mol2_file.h"
#include "pteros/core/pteros_error.h"

using namespace pteros;
using namespace std;

#ifndef USE_OPENBABEL

MOL2_file::MOL2_file(string &fname): VMD_molfile_plugin_wrapper(fname){
    plugin = molfile_plugins["mol2"];
}

void MOL2_file::do_write(const Selection &sel, const Mol_file_content &what){
    throw Pteros_error("In order to write MOL2 files you need to compile with OpenBabel support!");
}

#else

MOL2_file::MOL2_file(string &fname): Babel_wrapper(fname){ }

#endif



