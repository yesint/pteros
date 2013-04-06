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

#ifndef FORMAT_RECOGNITION_H_INCLUDED
#define FORMAT_RECOGNITION_H_INCLUDED
#include <string>

namespace pteros {

/// List of supported file formats identified by extension
enum FILE_FORMATS {PDB_FILE, GRO_FILE, TRR_FILE, XTC_FILE, DCD_FILE, PTTOP_FILE};

/// Takes file name and returns code of the file format.
/// Throws error if format is not recognized
FILE_FORMATS recognize_format(std::string& fname);

}

#endif // FORMAT_RECOGNITION_H_INCLUDED
