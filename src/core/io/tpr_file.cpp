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

#include "tpr_file.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

#ifdef USE_GROMACS
#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#endif

using namespace std;
using namespace pteros;
using namespace Eigen;


TPR_file::~TPR_file(){
}

#ifdef USE_GROMACS

void TPR_file::open(char open_mode)
{
    if(open_mode=='w')
        throw Pteros_error("TPR files could not be written!");
}

bool TPR_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){
    t_inputrec ir;    
    gmx_mtop_t mtop;
    t_topology top;
    t_state state;

    read_tpx_state(fname.c_str(), &ir, &state, &mtop);
    top = gmx_mtop_t_to_t_topology(&mtop);

    natoms = top.atoms.nr;

    if(frame){
        frame->coord.resize(natoms);
        if(what.coord()) frame->box.set_matrix(Map<Matrix3f>((float*)&state.box,3,3));
    }

    for(int i=0;i<natoms;++i){
        Atom at;
        int resi = top.atoms.atom[i].resind;
        at.name = *(top.atoms.atomname[i]);
        at.type_name = *(top.atoms.atomtype[i]);
        at.resid = top.atoms.resinfo[resi].nr;
        at.resname = *(top.atoms.resinfo[resi].name);
        at.mass = top.atoms.atom[i].m;
        at.charge = top.atoms.atom[i].q;
        char c = top.atoms.resinfo[resi].chainid;
        if(c!='\0') at.chain = c;
        if(top.atoms.pdbinfo){
            at.beta = top.atoms.pdbinfo[i].bfac;
            at.occupancy = top.atoms.pdbinfo[i].occup;
        }

        if(what.coord())
            for(int j=0;j<3;++j) frame->coord[i](j) = state.x[i][j];

        if(what.atoms()){
            // Add atoms to system
            append_atom_in_system(*sys,at);
        } else {
            // Update atoms in system
            atom_in_system(*sys,i) = at;
        }
    }

    return true;
}

#else

void TPR_file::open(char open_mode)
{
    throw Pteros_error("Pteros is compled without Gromacs support. Can't read TPR files!");
}


bool TPR_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){}

#endif
