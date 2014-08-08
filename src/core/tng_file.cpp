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

#include "tng_file.h"
#include "molfile_plugin.h"
#include "pteros/core/pteros_error.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

TNG_file::TNG_file(string fname, char open_mode): Mol_file(fname, open_mode){
    if(open_mode=='r'){
        if( tng_util_trajectory_open(fname.c_str(),'r',&trj) != TNG_SUCCESS )
            throw Pteros_error("Can't open TNG file '"+fname+"' for reading!");

        int64_t d_exp;
        tng_distance_unit_exponential_get(trj, &d_exp);
        unit_conv = std::pow(10,-9-d_exp);

        tng_num_particles_get(trj, &n_atoms);

        int64_t n1,n2;
        tng_num_frames_get(trj, &n1);
        tng_num_frame_sets_get(trj, &n2);
        cout << n1 << " " << n2 << endl;
        cur_fr = 1;

    } else {

    }
}

TNG_file::~TNG_file()
{
    tng_util_trajectory_close(&trj);
}

bool TNG_file::do_read(System *sys, Frame *frame, const Mol_file_content &what)
{
    bool stat;
    char datatype;
    int64_t len;
    double time_stamp;

    if(what.coordinates || what.trajectory){

        frame->coord.resize(n_atoms);

        float* ptr = (float*)&frame->coord.front();

        stat = tng_util_pos_read_range(trj,cur_fr,cur_fr, &ptr, &len);
        cout << cur_fr << " " << stat << " " << frame->coord[0].transpose() << endl;
        if(stat) cur_fr++;

    }

    return stat;
}

void TNG_file::do_write(const Selection &sel, const Mol_file_content &what)
{

}

