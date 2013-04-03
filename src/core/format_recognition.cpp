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

#include "pteros/core/format_recognition.h"
#include "pteros/core/pteros_error.h"

namespace pteros {

FILE_FORMATS recognize_format(std::string& fname){
    std::string ftype = fname.substr(fname.find_last_of(".") + 1);
    if(ftype=="xtc") return XTC_FILE;
    else if(ftype=="trr") return TRR_FILE;
    else if(ftype=="pdb") return PDB_FILE;
    else if(ftype=="gro") return GRO_FILE;
    else if(ftype=="top") return TOP_FILE;
    else if(ftype=="dcd") return DCD_FILE;
    else if(ftype=="tpr") return TPR_FILE;
    else if(ftype=="pttop") return PTTOP_FILE;
    else throw Pteros_error("File extension "+ftype+ " not recognized!");
}

}
