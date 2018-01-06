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


// Remove this define to use internal PDB reader.
/*
Internal reader is currently broken!
On write very large files fail to align coordinate records correctly
For now VMD molfile plugin is used for PDB files
*/

#define VMD_PDB

#ifndef PDB_FILE_H
#define PDB_FILE_H

#include "vmd_molfile_plugin_wrapper.h"

namespace pteros {

#ifdef VMD_PDB

/// Use VMD plugin for PDB
class PDB_file: public VMD_molfile_plugin_wrapper {
public:    

    PDB_file(std::string& fname);

    virtual Mol_file_content get_content_type() const {        
        return Mol_file_content().atoms(true).coord(true);
    }

};


#else
// Use internal PDB reader

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
This implementation is currently broken!
On write very large files fail to align coordinate records correctly
For now VMD molfile plugin is used for PDB files
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

class PDB_file: public Mol_file {
public:

    PDB_file(std::string fname, char open_mode);
    ~PDB_file();

    virtual bool do_read(System *sys, Frame *frame, Mol_file_content what);
    virtual void do_write(Selection &sel, Mol_file_content what);

    virtual Mol_file_content get_content_type(){
        Mol_file_content c;
        c.structure = true;
        c.coordinates = true;
        return c;
    }

protected:
    std::fstream f;    
};

#endif

}
#endif /* MOL_FILE_H */

