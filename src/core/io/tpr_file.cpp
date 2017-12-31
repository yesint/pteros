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

string coulomb_names[] = {"Cut-off", "Reaction-Field", "Generalized-Reaction-Field",
"PME", "Ewald", "P3M-AD", "Poisson", "Switch", "Shift", "User",
"Generalized-Born", "Reaction-Field-nec", "Encad-shift",
"PME-User", "PME-Switch", "PME-User-Switch",
"Reaction-Field-zero"};

string vdw_names[] = {"Cut-off", "Switch", "Shift", "User", "Encad-shift","PME"};

string mod_names[] = {"Potential-shift-Verlet", "Potential-shift", "None", "Potential-switch", "Exact-cutoff", "Force-switch"};


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
        at.type = top.atoms.atom[i].type;
        at.type_name = *(top.atoms.atomtype[i]);
        at.resid = top.atoms.resinfo[resi].nr;
        at.resname = *(top.atoms.resinfo[resi].name);
        at.mass = top.atoms.atom[i].m;
        at.charge = top.atoms.atom[i].q;
        at.element_number = top.atoms.atom[i].atomnumber;
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
        Force_field& ff = sys->get_force_field();
        ff.bonds.clear();

        ff.natoms = natoms;

        map<int,bool> is_bond_type;
        map<int,bool> is_lj14_type;
        map<int,bool> is_settle_type;
        map<int,int> lj14_type_map;

        vector<float> c6_l, c12_l;
        int atnr = top.idef.atnr;

        // lists for c6 and c12
        c6_l.reserve(atnr*atnr);
        c12_l.reserve(atnr*atnr);

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
            case F_LJ:
                c6_l.push_back(top.idef.iparams[i].lj.c6);
                c12_l.push_back(top.idef.iparams[i].lj.c12);
                break;
            case F_LJ14:
                is_lj14_type[i]=true;
                ff.LJ14_interactions.push_back({top.idef.iparams[i].lj14.c6A,top.idef.iparams[i].lj14.c12A});
                lj14_type_map[i] = ff.LJ14_interactions.size()-1;
                break;
            default:
                is_bond_type[i]=false;
                is_settle_type[i]=false;
                is_lj14_type[i]=false;
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
                        ff.bonds.push_back({a1,a2});
                    }
                }

                else if(is_settle_type[top.idef.il[it].iatoms[0]]){
                    for (int i=0; i<top.idef.il[it].nr; ){
                        type = top.idef.il[it].iatoms[i++];
                        a1 = top.idef.il[it].iatoms[i++];
                        a2 = top.idef.il[it].iatoms[i++];
                        a3 = top.idef.il[it].iatoms[i++];
                        // One settles entry is *two* O-H bonds!

                        ff.bonds.push_back({a1,a2});
                        ff.bonds.push_back({a1,a3});
                    }
                }

                else if(is_lj14_type[top.idef.il[it].iatoms[0]]){
                    int Nlj14 = ff.LJ14_interactions.size();
                    for (int i=0; i<top.idef.il[it].nr; ){
                        type = top.idef.il[it].iatoms[i++];
                        a1 = top.idef.il[it].iatoms[i++];
                        a2 = top.idef.il[it].iatoms[i++];
                        if(a1>a2) std::swap(a1,a2);
                        ff.LJ14_pairs[a1*natoms+a2]=lj14_type_map[type];
                    }
                }

            }
        }

        // Non-bond idefs
        // Here is how access is given in GROMACS between atoms with atomtypes A and B:
        // float c6=mtop->ffparams.iparams[B*atnr+A].lj.c6;
        // float c12=mtop->ffparams.iparams[B*atnr+A].lj.c12;

        ff.LJ_C6.resize(atnr,atnr);
        ff.LJ_C12.resize(atnr,atnr);

        for(int A=0; A<atnr; ++A){
            for(int B=0; B<atnr; ++B){
                ff.LJ_C6(A,B) = c6_l[B*atnr+A];
                ff.LJ_C6(B,A) = c6_l[B*atnr+A];
                ff.LJ_C12(A,B) = c12_l[B*atnr+A];
                ff.LJ_C12(B,A) = c12_l[B*atnr+A];
            }
        }


        // Molecules
        for(int i=0; i<top.mols.nr; ++i){
            ff.molecules.push_back({top.mols.index[i],top.mols.index[i+1]-1});
        }

        // Exclusions
        ff.exclusions.resize(natoms);
        int b,e;
        for(int i=0; i<top.excls.nr; ++i){
            //cout << top.excls.index[i] << " " << top.excls.index[i+1]-1 << endl;
            b = top.excls.a[top.excls.index[i]];
            e = top.excls.a[top.excls.index[i+1]-1];
            for(int j=b; j<e; ++j){
                //cout << i << " " << j << endl;
                if(j!=i) ff.exclusions[i].insert(j);
            }
            //LOG()->info("{}:{}",top.excls.a[top.excls.index[i]],top.excls.a[top.excls.index[i+1]-1]);
        }

        ff.fudgeQQ = top.idef.fudgeQQ;
        ff.rcoulomb = ir.rcoulomb;
        ff.epsilon_r = ir.epsilon_r;
        ff.epsilon_rf = ir.epsilon_rf;
        ff.rcoulomb_switch= ir.rcoulomb_switch;
        ff.rvdw_switch= ir.rvdw_switch;
        ff.rvdw= ir.rvdw;
        ff.coulomb_type= coulomb_names[ir.coulombtype];
        ff.coulomb_modifier= mod_names[ir.coulomb_modifier];
        ff.vdw_type= vdw_names[ir.vdwtype];
        ff.vdw_modifier= mod_names[ir.vdw_modifier];

        ff.setup_kernels(); // Setup kernels
        ff.ready = true;

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
