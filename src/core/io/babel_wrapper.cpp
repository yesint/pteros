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


#include "babel_wrapper.h"
#include "pteros/core/pteros_error.h"
#include "openbabel/atom.h"
#include "openbabel/residue.h"
#include "openbabel/bondtyper.h"
#include "openbabel/generic.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

using namespace std;
using namespace pteros;
using namespace Eigen;


Babel_wrapper::Babel_wrapper(string &fname): Mol_file(fname){ }

void Babel_wrapper::open(char open_mode)
{
    string ext = fname.substr(fname.find_last_of(".") + 1);
    bool ok;

    if(open_mode=='r'){
        ok = conv.SetInFormat(ext.c_str());
        if(!ok) throw Pteros_error("OpenBabel can't read format '{}'!",ext);
    } else {
        ok = conv.SetOutFormat(ext.c_str());
        if(!ok) throw Pteros_error("OpenBabel can't write format '{}'!",ext);
    }
}

Babel_wrapper::~Babel_wrapper()
{

}

bool Babel_wrapper::do_read(System *sys, Frame *frame, const Mol_file_content &what)
{
    bool ok;
    ok = conv.ReadFile(&mol,fname);
    if(!ok) throw Pteros_error("Babel can't read file {}",fname);

    int natoms = mol.NumAtoms();

    if(what.atoms()){
        // Allocate atoms in the system
        allocate_atoms_in_system(*sys, natoms);

        Atom at;
        // Cycle over atoms
        int i = 0;
        FOR_ATOMS_OF_MOL(ba, mol)
        {
            auto res = ba->GetResidue();

            at.name = res->GetAtomID(&(*ba));
            at.resname = res->GetName();
            at.resid = res->GetNum();
            at.chain = res->GetChain();
            at.occupancy = 0.0;
            at.beta = 0.0;
            at.charge = ba->GetPartialCharge();
            at.element_number = ba->GetAtomicNum();
            at.mass = ba->GetAtomicMass();

            set_atom_in_system(*sys,i,at);
            ++i;
        }
        sys->assign_resindex();
    }

    if(what.coord()){
        frame->coord.resize(natoms);
        int i=0;
        FOR_ATOMS_OF_MOL(ba, mol)
        {
            frame->coord[i](0) = ba->GetX() * 0.1;
            frame->coord[i](1) = ba->GetY() * 0.1;
            frame->coord[i](2) = ba->GetZ() * 0.1;
            ++i;
        }
    }

    return true;
}

void Babel_wrapper::do_write(const Selection &sel, const Mol_file_content &what)
{
    if(what.atoms() && what.coord()){
        selection_to_obmol(sel,mol);
    }

    conv.WriteFile(&mol,fname);
}
