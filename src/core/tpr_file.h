/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
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

#ifndef TPR_FILE_H
#define TPR_FILE_H

#include <string>
#include <fstream>
#include "pteros/core/mol_file.h"
#include "pteros/core/pteros_error.h"

namespace pteros {


class TPR_file: public Mol_file {
public:

    TPR_file(std::string fname, char mode);
    ~TPR_file();   

    virtual Mol_file_content get_content_type(){
        Mol_file_content c;
        c.topology = true;
        return c;
    }

protected:
    std::string dump_name;
    char open_mode;
    std::string file_name;

    virtual void do_write(Selection &sel, Mol_file_content what){
        throw Pteros_error("TPR files could not be written by Pteros! Use GROMACS to make them.");
    }

    virtual bool do_read(System *sys, Frame *frame, Mol_file_content what);
};

}
#endif /* MOL_FILE_H */
