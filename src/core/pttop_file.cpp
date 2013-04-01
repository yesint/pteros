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

#include "pttop_file.h"
#include "pteros/core/pteros_error.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


PTTOP_file::PTTOP_file(string fname, char mode): Mol_file(fname, mode)
{    
    file_name = fname;
    open_mode = mode;
}

PTTOP_file::~PTTOP_file(){
}

bool PTTOP_file::do_read(System *sys, Frame *frame, Mol_file_content what){
    if(open_mode=='r'){

    } else {
        throw Pteros_error("PTTOP files could not be written from Pteros!");
    }

}
