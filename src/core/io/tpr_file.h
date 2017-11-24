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

#pragma once

#include "pteros/core/mol_file.h"

namespace pteros {


class TPR_file: public Mol_file {
public:
    TPR_file(std::string& fname): Mol_file(fname) {}
    virtual void open(char open_mode);
    ~TPR_file();

    virtual Mol_file_content get_content_type() const {        
        return Mol_file_content().atoms(true).coord(true).top(true);
    }

protected:        

    virtual void do_write(const Selection &sel, const Mol_file_content& what) {}

    virtual bool do_read(System *sys, Frame *frame, const Mol_file_content& what);
};

}

