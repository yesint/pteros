/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2014, Semen Yesylevskyy
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

#ifndef PTTOP_FILE_H
#define PTTOP_FILE_H

#include <string>
#include <fstream>
#include "pteros/core/mol_file.h"
#include "pteros/core/pteros_error.h"

namespace pteros {


class PTTOP_file: public Mol_file {
public:

    PTTOP_file(std::string fname, char mode);
    ~PTTOP_file();

    virtual Mol_file_content get_content_type(){
        Mol_file_content c;
        c.topology = true;
        c.structure= true;
        c.coordinates= true;
        return c;
    }

protected:    

    virtual void do_write(Selection &sel, Mol_file_content what){
        throw Pteros_error("PTTOP files could be produced by the dedicated script only!");
    }

    virtual bool do_read(System *sys, Frame *frame, Mol_file_content what);

    std::string file_name;
    char open_mode;
};

}
#endif /* MOL_FILE_H */
