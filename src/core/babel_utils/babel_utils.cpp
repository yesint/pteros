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

        oba->SetAtomicNum(at.element_number);
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

    mol.SetPartialChargesPerceived(); // Needed again (openbabel is written by retarded perverts!!!)
}
