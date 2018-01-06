/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *  
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#include "xtc_file.h"

using namespace std;
using namespace pteros;

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

