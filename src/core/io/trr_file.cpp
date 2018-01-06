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


#include "trr_file.h"

using namespace std;
using namespace pteros;


int TRR_file::read_num_atoms(char* fname, int* num){
    return read_trr_natoms(fname, num);
}

int TRR_file::read_record(XDRFILE *xd, int natoms, int *step,
                        float *time, matrix box,rvec *x){
    float lambda;
    return read_trr(xd,natoms,step,time,&lambda,box,x,NULL,NULL);
}

int TRR_file::write_record(XDRFILE *xd, int natoms, int step,
                           float time, matrix box, rvec *x){
    float lambda = 1.0;
    return write_trr(xd,natoms,step,time,lambda,box,x,NULL,NULL);
}

