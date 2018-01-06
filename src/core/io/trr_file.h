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


#ifndef TRR_FILE_H
#define TRR_FILE_H

#include "gromacs_trajectory_file.h"

namespace pteros {

/** TRR reader
  */
class TRR_file: public Gromacs_trajectory_file {
public:
    TRR_file(std::string fname): Gromacs_trajectory_file(fname) {}

protected:

    virtual int read_num_atoms(char* fname, int* num);
    virtual int read_record(XDRFILE *xd, int natoms, int *step,
                            float *time, matrix box,rvec *x);
    virtual int write_record(XDRFILE *xd, int natoms, int step,
                             float time, matrix box, rvec *x);
};

}
#endif // GROMACS_TRAJECTORY_H

