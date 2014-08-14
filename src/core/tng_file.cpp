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

using namespace std;
using namespace pteros;
using namespace Eigen;

extern molfile_plugin_t tng_plugin;

TNG_file::TNG_file(string fname, char open_mode): VMD_molfile_plugin_wrapper(fname,open_mode){
    plugin = &tng_plugin;
    accepted_format = TNG_FILE;
    open(fname,open_mode);
}

/*
TNG_file::TNG_file(string fname, char open_mode): Mol_file(fname, open_mode){
    if(open_mode=='r'){
        if( tng_util_trajectory_open(fname.c_str(),'r',&trj) != TNG_SUCCESS )
            throw Pteros_error("Can't open TNG file '"+fname+"' for reading!");

        int64_t d_exp;
        tng_distance_unit_exponential_get(trj, &d_exp);
        unit_conv = std::pow(10,-9-d_exp);

        tng_num_particles_get(trj, &n_atoms);

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

    void* values = 0;
    int64_t frame_num;
    double frame_time;

    if(what.coordinates || what.trajectory){

        stat = tng_util_particle_data_next_frame_read(trj, TNG_TRAJ_POSITIONS, &values,
                                                      &datatype, &frame_num, &frame_time);
        cout << cur_fr << " " << stat << endl;
        if(stat){
            cur_fr++;
        }
    }

    return stat;
}

void TNG_file::do_write(const Selection &sel, const Mol_file_content &what)
{

}

*/
