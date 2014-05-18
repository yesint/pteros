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

#include "gromacs_trajectory_file.h"
#include <iostream>
#include "pteros/core/pteros_error.h"
#include <stdlib.h>

using namespace std;
using namespace pteros;

Gromacs_trajectory_file::Gromacs_trajectory_file(std::string fname, char openmode):
    Mol_file(fname,openmode)
{    
}

void Gromacs_trajectory_file::open(std::string fname, char openmode){
    mode = openmode;
    if(mode=='r')
        cout << "Opening trajectory file " << fname << " for reading..." << endl;
    else
        cout << "Opening trajectory file " << fname << " for writing..." << endl;

    xd = xdrfile_open(fname.c_str(),&mode);
    if(!xd) throw Pteros_error("Can't open trajectory file "+fname+" with mode "+mode);

    if(mode == 'r'){ //READING
        // Get number of atoms
        ret = read_num_atoms(const_cast<char*>(fname.c_str()), &natoms);

        if(ret == exdrOK)
            cout << "There are " << natoms << " atoms in this trajectory" << endl;
        else
            throw Pteros_error("Failed to get number of atoms in trajectory file "+fname);

        // Allocate data for trajectory frame
        //x = (rvec*) calloc(natoms,sizeof(*x));
        //cout << "Memory buffer allocated"<<endl;

    } else {
        // For writing number of atoms will be supplied later on first write operation
        // Data will also be allocated there.

        // Prepare the box just in case
        box[0][0] = 0.0; box[0][1] = 0.0; box[0][2] = 0.0;
        box[1][0] = 0.0; box[1][1] = 0.0; box[1][2] = 0.0;
        box[2][0] = 0.0; box[2][1] = 0.0; box[2][2] = 0.0;
    }

    fr = 0; //First frame
    x = NULL;
    time_per_frame = 0;
}

Gromacs_trajectory_file::~Gromacs_trajectory_file(){
    if(x){
        free(x);
        x = NULL;
    }

    if(xd){
        xdrfile_close(xd);
        xd = NULL;
    }

    if(mode=='r')
        cout << "Raw number of frames read from file: " << fr
             << ". Each frame takes " << time_per_frame/(float)CLOCKS_PER_SEC/(float)fr << " s."
             <<  " Total reading time " << time_per_frame/(float)CLOCKS_PER_SEC << " s." <<endl;
    else
        cout << "Raw number of frames written to file: " << fr
             << ". Each frame takes " << time_per_frame/(float)CLOCKS_PER_SEC/(float)fr << " s."
             <<  " Total writing time " << time_per_frame/(float)CLOCKS_PER_SEC << " s." <<endl;
}

bool Gromacs_trajectory_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){
    clock_t t1 = clock();
    int i;
    float lambda;    

    frame->coord.resize(natoms);
    ret = read_record(xd,natoms,&step,&frame->time, box, (rvec*)&frame->coord.front());

    if (ret != exdrENDOFFILE){
        if(ret != exdrOK){
            cout << "Trajectory seems to be is corrupted or incomplete" << endl;
            return false;
        } else {
            // Everything is Ok
            ++fr;
            // Get box
            Eigen::Matrix3f b;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    b(i,j) = box[i][j];
            frame->box.modify(b);
        }
    } else {
        return false; // End of file
    }


    time_per_frame += clock()-t1;    
    return true; // Allows to proceed to next frame
}

void Gromacs_trajectory_file::do_write(const Selection &sel, const Mol_file_content &what){
    clock_t t1 = clock();

    if(!x){
        // Do initialization on first step
        natoms = sel.size();        
        // Allocate data for trajectory frame
        x = (rvec*) calloc(natoms,sizeof(*x));
        fr = 0;
    }

    // Copy data to internal storage
    //cout << "----------------------------"<<endl;
    for(int i=0;i<natoms;++i){
        x[i][0] = sel.X(i);
        x[i][1] = sel.Y(i);
        x[i][2] = sel.Z(i);

        //cout << i << ": " << sel.XYZ(i).transpose() << endl;
    }

    // Set box
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            box[i][j] = sel.get_system()->Box(sel.get_frame()).get_box()(i,j);


    // Write
    ret = write_record(xd,natoms,fr,sel.get_system()->Time(sel.get_frame()),box,x);    

    if(ret!=0) throw Pteros_error("Unable to write frame!");
    ++fr;

    time_per_frame += clock()-t1;
}
