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

#include "tpr_file.h"
#include "pteros/core/pteros_error.h"
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

using namespace std;
using namespace pteros;
using namespace Eigen;


TPR_file::TPR_file(string fname, char mode): Mol_file(fname, mode)
{
    dump_name = "_dump_of_tpr_for_pteros";
    open_mode = mode;
    file_name = fname;
}

TPR_file::~TPR_file(){
    remove(dump_name.c_str());
}

bool TPR_file::do_read(System *sys, Frame *frame, Mol_file_content what){
    if(open_mode=='r'){

        cout << "Running gmxdump to dump content of TPR file " << file_name << "..." << endl;
        string cmd("gmxdump -sys -s ");
        cmd += file_name + " > " + dump_name;
        int ret = system(cmd.c_str());
        if(ret>0) throw Pteros_error("Error executing gmxdump! Check output above for error message.");
    } else {
        throw Pteros_error("TPR files could not be written from Pteros!");
    }


}
