/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2023, Semen Yesylevskyy
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


#include "selection_to_obmol.h"
#include "openbabel/atom.h"
#include "openbabel/residue.h"
#include "openbabel/bondtyper.h"
#include "openbabel/generic.h"
#include "openbabel/obfunctions.h"
#include "openbabel/obiter.h"

using namespace std;
using namespace pteros;

namespace pteros {
// Converts Selection to OBMol
void selection_to_obmol(const Selection& sel, OpenBabel::OBMol &mol, bool babel_bonds)
{
    mol.Clear();

    // map of residues
    map<int,OpenBabel::OBResidue*> reslist;

    mol.BeginModify();
    mol.SetPartialChargesPerceived(); // Needed to suppress charges computations

    for(int i=0;i<sel.size();++i){
        auto& at = sel.atom(i);

        // Create new atom in this mol
        auto oba = mol.NewAtom();

        oba->SetAtomicNum(at.atomic_number);
        oba->SetPartialCharge(at.charge);
        oba->SetVector(10.0*sel.x(i),10.0*sel.y(i),10.0*sel.z(i));

        // Create new residue if needed
        if(reslist.count(at.resid)==0){
            OpenBabel::OBResidue* obr = mol.NewResidue();
            obr->SetNum(at.resid);
            obr->SetChain(at.chain);
            obr->SetName(at.resname);
            reslist[at.resid] = obr;
        }

        reslist[at.resid]->AddAtom(oba);
        reslist[at.resid]->SetAtomID(oba,at.name);
    }

    if(babel_bonds){
        mol.ConnectTheDots();
        // Guess bond orders and aromaticity
        mol.PerceiveBondOrders();
    } else {
        // Get bonds from pteros
        auto con = sel.get_internal_bonds(0.18,false); // Non-periodic by default
        // Set bonds manually
        for(int i=0; i<con.size(); ++i){
            for(int j=0; j<con[i].size(); ++j) mol.AddBond(i,con[i][j],1);
        }
    }

    FOR_ATOMS_OF_MOL(matom, mol)
        OBAtomAssignTypicalImplicitHydrogens(&*matom);

    mol.EndModify();

    // Needed again (openbabel is written by retarded perverts!!!)
    mol.SetPartialChargesPerceived();
}

}





