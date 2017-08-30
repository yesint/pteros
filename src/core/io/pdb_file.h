/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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
