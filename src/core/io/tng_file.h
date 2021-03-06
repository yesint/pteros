/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2020, Semen Yesylevskyy
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



#ifndef TNG_FILE_H
#define TNG_FILE_H

#include <string>
#include "vmd_molfile_plugin_wrapper.h"
#include "tng/tng_io.h"

namespace pteros {

/// Use VMD plugin for TNG
class TNG_file: public VMD_molfile_plugin_wrapper {
public:
    TNG_file(std::string& fname);

    virtual Mol_file_content get_content_type() const {                
        return Mol_file_content().atoms(true).traj(true);
    }

};

}
#endif /* MOL_FILE_H */


