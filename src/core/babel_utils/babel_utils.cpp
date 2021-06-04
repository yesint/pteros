/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2021, Semen Yesylevskyy
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


#include "babel_utils.h"
#include <openbabel/obfunctions.h>
#include <openbabel/obiter.h>

using namespace pteros;
using namespace std;

void pteros::selection_to_obmol(const Selection& sel, OpenBabel::OBMol &mol, bool babel_bonds)
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
