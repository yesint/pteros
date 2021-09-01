/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
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

#include "tpr_file.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "system_builder.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/idef.h"

// Generated version infor file
#include "gromacs_version_info.h"

#if (GROMACS_VERSION > 2020)
#include "gromacs/mdtypes/md_enums.h"
#endif


using namespace std;
using namespace pteros;
using namespace Eigen;


string coulomb_names[] = {"Cut-off","Reaction-Field", "Generalized-Reaction-Field",
                          "PME", "Ewald", "P3M-AD", "Poisson", "Switch", "Shift", "User",
                          "Generalized-Born", "Reaction-Field-nec", "Encad-shift",
                          "PME-User", "PME-Switch", "PME-User-Switch",
                          "Reaction-Field-zero"};

string vdw_names[] = {"Cut-off", "Switch", "Shift", "User", "Encad-shift","PME"};

string mod_names[] = {"Potential-shift-Verlet", "Potential-shift", "None", "Potential-switch",
                      "Exact-cutoff", "Force-switch"};


void TprFile::open(char open_mode)
{
    if(open_mode=='w')
        throw PterosError("TPR files could not be written!");
}

void TprFile::close()
{

}

bool TprFile::do_read(System *sys, Frame *frame, const FileContent &what){
    t_inputrec ir;    
    gmx_mtop_t mtop;
    t_topology top;
    t_state state;

    read_tpx_state(fname.c_str(), &ir, &state, &mtop);

    // Only works for Gromacs>=2018.x
    top = gmx_mtop_t_to_t_topology(&mtop,false);

    natoms = top.atoms.nr;

    if(frame){
        frame->coord.resize(natoms);
        if(what.coord()) frame->box.set_matrix(Map<Matrix3f>((float*)&state.box,3,3));
    }

    SystemBuilder builder(sys);

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
        at.atomic_number = top.atoms.atom[i].atomnumber;
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
            builder.add_atom(at);
        } else {
            // Update atoms in system
            builder.atom(i) = at;
        }
    }

    if(what.atoms()) sys->assign_resindex();

    // Read topology
    if(what.top()){
        ForceField& ff = sys->get_force_field();
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
                ff.LJ14_interactions.emplace_back(top.idef.iparams[i].lj14.c6A,top.idef.iparams[i].lj14.c12A);
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
                        ff.bonds.emplace_back(a1,a2);
                    }
                }

                else if(is_settle_type[top.idef.il[it].iatoms[0]]){
                    for (int i=0; i<top.idef.il[it].nr; ){
                        type = top.idef.il[it].iatoms[i++];
                        a1 = top.idef.il[it].iatoms[i++];
                        a2 = top.idef.il[it].iatoms[i++];
                        a3 = top.idef.il[it].iatoms[i++];
                        // One settles entry is *two* O-H bonds!
                        ff.bonds.emplace_back(a1,a2);
                        ff.bonds.emplace_back(a1,a3);
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
            ff.molecules.emplace_back(top.mols.index[i],top.mols.index[i+1]-1);
        }

        // Exclusions

        // Starting from Gromacs 2021 exclusions are removed from global top
#if (GROMACS_VERSION <= 2020)
        ff.exclusions.resize(natoms);
        int b,e;
        for(int i=0; i<top.excls.nr; ++i){
            b = top.excls.a[top.excls.index[i]];
            e = top.excls.a[top.excls.index[i+1]-1];
            for(int j=b; j<e; ++j){
                //cout << i << " " << j << endl;
                if(j!=i) ff.exclusions[i].insert(j);
            }
        }
#else

        // Extract exclusions from molecular blocks
        ff.exclusions.resize(natoms);
        // Running global index of atoms
        size_t global_counter = 0;
        // Cycle over molecular blocks
        for(size_t bl=0; bl<mtop.molblock.size();++bl){
            // Mol type
            size_t mol_t = mtop.molblock[bl].type;
            // Cycle over molecules in this block
            for(size_t mol=0; mol<mtop.molblock[bl].nmol; ++mol){
                // cycle over atoms in one molecule
                for(size_t a=0; a<mtop.moltype[mol_t].excls.size(); ++a){
                    // Go over the list of excluded local indexes for this atom
                    for(size_t ind=0; ind<mtop.moltype[mol_t].excls[a].size(); ++ind){
                        // Global excluded indes for local atom a
                        size_t global_excl = mtop.moltype[mol_t].excls[a][ind] + global_counter;
                        // Add this exclusion
                        ff.exclusions[global_counter].insert(global_excl);
                    }
                    // Increment global atom counter
                    ++global_counter;
                }
            }
        }

#endif

        ff.fudgeQQ = top.idef.fudgeQQ;
        ff.rcoulomb = ir.rcoulomb;
        ff.epsilon_r = ir.epsilon_r;
        ff.epsilon_rf = ir.epsilon_rf;
        ff.rcoulomb_switch= ir.rcoulomb_switch;
        ff.rvdw_switch= ir.rvdw_switch;
        ff.rvdw= ir.rvdw;

        // Since Gromacs 2021 interaction types are named enums, so cast to int
        ff.coulomb_type= coulomb_names[int(ir.coulombtype)];
        ff.coulomb_modifier= mod_names[int(ir.coulomb_modifier)];
        ff.vdw_type= vdw_names[int(ir.vdwtype)];
        ff.vdw_modifier= mod_names[int(ir.vdw_modifier)];

        ff.setup_kernels(); // Setup kernels
        ff.ready = true;

    } // topology

    return true;
}



