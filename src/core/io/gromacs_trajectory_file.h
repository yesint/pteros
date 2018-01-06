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

