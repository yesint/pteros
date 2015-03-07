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

#ifndef XTC_FILE_H
#define XTC_FILE_H

#include "gromacs_trajectory_file.h"

namespace pteros {

/** TRR reader
  */
class XTC_file: public Gromacs_trajectory_file {
public:
    XTC_file(std::string fname): Gromacs_trajectory_file(fname) {}

protected:
    virtual int read_num_atoms(char* fname, int* num);
    virtual int read_record(XDRFILE *xd, int natoms, int *step,
                            float *time, matrix box,rvec *x);
    virtual int write_record(XDRFILE *xd, int natoms, int step,
                             float time, matrix box, rvec *x);
};

}
#endif // GROMACS_TRAJECTORY_H
