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
#include "gromacs/topology/idef.h"
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

    // Read atoms and coordinates
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

    // Read topology
    if(what.top()){
        ff_in_system(*sys).bonds.clear();

        ff_in_system(*sys).fudgeQQ = top.idef.fudgeQQ;

        map<int,bool> is_bond_type;
        map<int,bool> is_settle_type;

        //LOG()->info("ntypes={}, F_NRE={}",top.idef.ntypes,F_NRE);
        for (int i = 0; i < top.idef.ntypes; i++) {
            switch(top.idef.functype[i]){
            case F_BONDS:
            case F_G96BONDS:
            case F_HARMONIC:
            case F_CONSTR:
            case F_CONSTRNC:
                is_bond_type[i]=true;
                break;
            case F_SETTLE:
                is_settle_type[i]=true;
                break;
            default:
                is_bond_type[i]=false;
                is_settle_type[i]=false;
                break;
            }
        }

        int type,a1,a2,a3;
        float r0,k;
        for(int it=0; it<F_NRE; ++it){ // Over all interaction terms
            if(top.idef.il[it].nr>0){

                if(is_bond_type[top.idef.il[it].iatoms[0]]){
                    for (int i=0; i<top.idef.il[it].nr; ){
                        type = top.idef.il[it].iatoms[i++];
                        a1 = top.idef.il[it].iatoms[i++];
                        a2 = top.idef.il[it].iatoms[i++];
                        //ff_in_system(*sys).bonds_per_atom[a1].push_back(a2);
                        //ff_in_system(*sys).bonds_per_atom[a2].push_back(a1);
                        ff_in_system(*sys).bonds.push_back({a1,a2});

                        //r0 = top.idef.iparams[type].harmonic.rA;
                        //k = top.idef.iparams[type].harmonic.krA;
                        //LOG()->info("a1={}, a2={}, r0={}, k={}",a1,a2,r0,k);
                        //LOG()->info("a1={}, a2={}, n1={}, n2={}",a1,a2,*(top.atoms.atomname[a1]),*(top.atoms.atomname[a2]));
                    }
                }

                else if(is_settle_type[top.idef.il[it].iatoms[0]]){
                    for (int i=0; i<top.idef.il[it].nr; ){
                        type = top.idef.il[it].iatoms[i++];
                        a1 = top.idef.il[it].iatoms[i++];
                        a2 = top.idef.il[it].iatoms[i++];
                        a3 = top.idef.il[it].iatoms[i++];
                        // One settles entry is *two* O-H bonds!
                        //ff_in_system(*sys).bonds_per_atom[a1].push_back(a2);
                        //ff_in_system(*sys).bonds_per_atom[a2].push_back(a1);

                        //ff_in_system(*sys).bonds_per_atom[a1].push_back(a3);
                        //ff_in_system(*sys).bonds_per_atom[a3].push_back(a1);
                        ff_in_system(*sys).bonds.push_back({a1,a2});
                        ff_in_system(*sys).bonds.push_back({a1,a3});


                        //LOG()->info("a1={}, a2={}, n1={}, n2={}",a1,a2,*(top.atoms.atomname[a1]),*(top.atoms.atomname[a2]));
                        //LOG()->info("a1={}, a2={}, n1={}, n2={}",a1,a3,*(top.atoms.atomname[a1]),*(top.atoms.atomname[a3]));
                    }
                }
            }
        }



    } // topology

    return true;
}

#else

void TPR_file::open(char open_mode)
{
    throw Pteros_error("Pteros is compled without Gromacs support. Can't read TPR files!");
}


bool TPR_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){}

#endif
