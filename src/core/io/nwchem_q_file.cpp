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

#include "nwchem_q_file.h"
#include "pteros/core/pteros_error.h"
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace pteros;
using namespace Eigen;

void Q_file::open(char open_mode)
{
    if(open_mode=='r'){
        f.open(fname.c_str(),ios_base::in);
        if(!f) throw Pteros_error("Can't open NWCHEM Q file '"+fname+"'' for reading");
    } else {
        f.open(fname.c_str(),ios_base::out);
        if(!f) throw Pteros_error("Can't open NWCHEM Q file '"+fname+"'' for writing");
    }
}

Q_file::~Q_file(){
    if(f){
        f.close();
    }
}

bool Q_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){

    int N,i;
    // tmp atom
    Atom tmp_atom;
    // Temporary coordinates
    Vector3f tmp_coor;

    // Read number of atoms
    f >> N >> i; // i is not used

    frame->coord.resize(N);

    tmp_atom.resname = "X";
    tmp_atom.chain = 'X';
    tmp_atom.beta = 0.0;
    tmp_atom.occupancy = 0.0;
    tmp_atom.type = -1;

    // Read coordinates
    for(i=0;i<N;++i){
        f >> tmp_atom.name >> tmp_coor(0) >> tmp_coor(1) >> tmp_coor(2) >> tmp_atom.charge;

        // Coordinates are in Angstrom, so we need to convert
        //tmp_coor /= 10.0;

        if(what & MFC_ATOMS){
            // Assign mass
            tmp_atom.mass = get_mass_from_atom_name(tmp_atom.name);
            // Add new atom to the system
            append_atom_in_system(*sys,tmp_atom);
        }

        if(what & MFC_COORD){
            // Add column of coordinates
            frame->coord[i] = tmp_coor;
        }
    }

    return true;
}

void Q_file::do_write(const Selection &sel, const Mol_file_content &what){
    throw Pteros_error("NWCHEM Q files are read-only!");
}
