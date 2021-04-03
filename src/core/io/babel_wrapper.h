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
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/file_handler.h"

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

namespace pteros {

/// Generic API for reading and writing any molecule file formats
class BabelWrapper: public FileHandler {
public:
    // High-level API        
    BabelWrapper(std::string& fname);
    virtual void open(char open_mode);
    virtual void close();

protected:       
    OpenBabel::OBConversion conv;
    OpenBabel::OBMol mol;

    // Tells if the format need bonds to be present
    virtual bool need_bonds() = 0;

    virtual bool do_read(System *sys, Frame *frame, const FileContent& what);
    virtual void do_write(const Selection &sel, const FileContent& what);
};

}





