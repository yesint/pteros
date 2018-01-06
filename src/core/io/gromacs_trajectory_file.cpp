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


#include "gromacs_trajectory_file.h"
#include <iostream>
#include "pteros/core/pteros_error.h"
#include <stdlib.h>
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;

Gromacs_trajectory_file::Gromacs_trajectory_file(string &fname): Mol_file(fname)
{
    x = nullptr;
    xd = nullptr;
}

void Gromacs_trajectory_file::open(char openmode){
    mode = openmode;

    xd = xdrfile_open(fname.c_str(),&mode);
    if(!xd) throw Pteros_error("Can't open trajectory file '{}' in mode '{}'",fname,mode);

    if(mode == 'r'){ //READING
        // Get number of atoms
        ret = read_num_atoms(const_cast<char*>(fname.c_str()), &natoms);

        if(ret != exdrOK)
            throw Pteros_error("Failed to get number of atoms in trajectory file {}",fname);

    } else {
        // For writing number of atoms will be supplied later on first write operation
        // Data will also be allocated there.

        // Prepare the box just in case
        box[0][0] = 0.0; box[0][1] = 0.0; box[0][2] = 0.0;
        box[1][0] = 0.0; box[1][1] = 0.0; box[1][2] = 0.0;
        box[2][0] = 0.0; box[2][1] = 0.0; box[2][2] = 0.0;
    }

    fr = 0; //First frame
    x = nullptr;
}

Gromacs_trajectory_file::~Gromacs_trajectory_file(){
    if(x){
        free(x);
        x = nullptr;
    }

    if(xd){
        xdrfile_close(xd);
        xd = nullptr;
    }
}

bool Gromacs_trajectory_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){
    int i;    

    frame->coord.resize(natoms);
    ret = read_record(xd,natoms,&step,&frame->time, box, (rvec*)&frame->coord.front());

    if (ret != exdrENDOFFILE){
        if(ret != exdrOK){
            LOG()->error("Trajectory seems to be is corrupted or incomplete");
            return false;
        } else {
            // Everything is Ok
            ++fr;
            // Get box
            Eigen::Matrix3f b;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    b(i,j) = box[j][i];
            frame->box.set_matrix(b);
        }
    } else {
        return false; // End of file
    }

    return true; // Allows to proceed to next frame
}

void Gromacs_trajectory_file::do_write(const Selection &sel, const Mol_file_content &what){
    if(!x){
        // Do initialization on first step
        natoms = sel.size();        
        // Allocate data for trajectory frame
        x = (rvec*) calloc(natoms,sizeof(*x));
        fr = 0;
    }

    // Copy data to internal storage    
    for(int i=0;i<natoms;++i){
        x[i][0] = sel.X(i);
        x[i][1] = sel.Y(i);
        x[i][2] = sel.Z(i);
    }

    // Set box
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            box[i][j] = sel.Box().get_matrix()(j,i);


    // Write
    ret = write_record(xd,natoms,fr,sel.get_system()->Time(sel.get_frame()),box,x);    

    if(ret!=0) throw Pteros_error("Unable to write frame!");
    ++fr;
}

