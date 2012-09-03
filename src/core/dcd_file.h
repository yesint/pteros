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

#ifndef DCD_FILE_H
#define DCD_FILE_H

#include <string>
#include "vmd_molfile_plugin_wrapper.h"

namespace pteros {

/// Generic API for reading and writing any molecule file formats
class DCD_file: public VMD_molfile_plugin_wrapper {
public:
    virtual Mol_file_content get_content_type(){
        Mol_file_content c;
        c.trajectory = true;
        return c;
    }

    DCD_file(std::string fname, char open_mode);
};

}
#endif /* MOL_FILE_H */
