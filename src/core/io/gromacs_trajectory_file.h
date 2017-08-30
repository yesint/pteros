/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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

#ifndef GROMACS_TRAJECTORY_FILE_H
#define GROMACS_TRAJECTORY_FILE_H

#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_trr.h"

#include "pteros/core/mol_file.h"

namespace pteros {

/** Base class for XTC and TRR readers
  */
class Gromacs_trajectory_file: public Mol_file {
    public:      

        Gromacs_trajectory_file(std::string& fname);

        virtual void open(char openmode);

        virtual ~Gromacs_trajectory_file();        

        virtual Mol_file_content get_content_type() const {
            return Mol_file_content().traj(true);
        }

    protected:

        XDRFILE* xd;
        int step, ret;
        float prec;
        matrix box;
        rvec* x;
        char mode;
        int fr;        

        virtual int read_num_atoms(char* fname, int* num) = 0;
        virtual int read_record(XDRFILE *xd, int natoms, int *step,
                                float *time, matrix box,rvec *x) = 0;
        virtual int write_record(XDRFILE *xd,
                      int natoms,int step,float time,
                      matrix box,rvec *x) = 0;

        /// Reads next frame into internal storage. Returns false if failed or EOF reached.
        virtual bool do_read(System *sys, Frame *frame, const Mol_file_content& what);
        /// Writes data from given frame to trajectory
        virtual void do_write(const Selection &sel, const Mol_file_content& what);
};

}
#endif // GROMACS_TRAJECTORY_H
