/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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

#include "xtc_file.h"

using namespace std;
using namespace pteros;

XTC_file::XTC_file(std::string fname, char openmode):
    Gromacs_trajectory_file(fname,openmode)
{
    open(fname,openmode);
}

int XTC_file::read_num_atoms(char* fname, int* num){
    return read_xtc_natoms(fname, num);
}

int XTC_file::read_record(XDRFILE *xd, int natoms, int *step,
                        float *time, matrix box,rvec *x){
    return read_xtc(xd,natoms,step,time,box,x,&prec);
}

int XTC_file::write_record(XDRFILE *xd, int natoms, int step,
                           float time, matrix box, rvec *x){
    int prec = 1000;
    return write_xtc(xd,natoms,step,time,box,x,prec);
}
