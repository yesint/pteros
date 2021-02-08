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


#include "pteros/core/mol_file.h"
#include "pteros/core/pteros_error.h"

#include "pdb_file.h"
#include "dcd_file.h"
#include "gro_file.h"
#include "mol2_file.h"
#include "xyz_file.h"
#include "trr_file.h"
#include "xtc_file.h"

#ifdef USE_TNGIO
#include "tng_file.h"
#endif

#ifdef USE_OPENBABEL
#include "pdbqt_file.h"
#endif

#ifdef USE_GROMACS
#include "tpr_file.h"
#endif

using namespace std;
using namespace pteros;

Mol_file::Mol_file(string& file_name){    
    fname = file_name;
}

Mol_file::~Mol_file(){    
}

bool Mol_file::read(System *sys, Frame *frame, const Mol_file_content &what){
    sanity_check_read(sys,frame,what);    
    return do_read(sys,frame,what);
}

void Mol_file::write(const Selection &sel, const Mol_file_content &what) {
    sanity_check_write(sel,what);
    do_write(sel,what);
}

void Mol_file::sanity_check_read(System *sys, Frame *frame, const Mol_file_content& what) const {
    auto c = get_content_type();
    if( !c.atoms() && what.atoms() )
        throw Pteros_error("Can't read structure from this file type!");
    if( !c.coord() && what.coord() )
        throw Pteros_error("Can't read coordinates from this file type!");
    if( !c.traj() && what.traj() )
        throw Pteros_error("Can't read coordinates from this file type!");
    if( !c.top() && what.top() )
        throw Pteros_error("Can't read topology from this file type!");
    if( !what.atoms() && !what.coord() && !what.traj() && !what.top() )
        throw Pteros_error("Nothing to read!");
    if(what.atoms() && !sys)
        throw Pteros_error("System should be provided to read structure!");
    if(what.top() && !sys)
        throw Pteros_error("System should be provided to read topology!");
    if((what.coord() || what.traj()) && !frame)
        throw Pteros_error("Frame should be provided to read coordinates!");
    if(what.atoms() && sys && sys->num_atoms()>0)
        throw Pteros_error("Can't read structure because system already has atoms!");    
}

void Mol_file::sanity_check_write(const Selection &sel, const Mol_file_content &what) const{
    if(!what.atoms() && !what.coord() && !what.traj() && !what.top())
        throw Pteros_error("Nothing to write!");
    if( !get_content_type().atoms() && what.atoms() )
        throw Pteros_error("Can't write structure from this file type!");
    if( !get_content_type().coord() && what.coord() )
        throw Pteros_error("Can't write coordinates from this file type!");
    if( !get_content_type().traj() && what.traj() )
        throw Pteros_error("Can't write coordinates from this file type!");
    if( !get_content_type().top() && what.top() )
        throw Pteros_error("Can't write topology from this file type!");
}

namespace pteros {

void get_element_from_atom_name(string& name, int &anum, float &mass){
    // Find first character, which is not digit to account for cases like 21C2
    int i = name.find_first_not_of("1234567890");

         if(name[i]=='C') { mass = 12.0; anum = 6; }
    else if(name[i]=='O') { mass = 16.0; anum = 8; }
    else if(name[i]=='N') { mass = 14.0; anum = 7; }
    else if(name[i]=='S') { mass = 32.0; anum = 16; }
    else if(name[i]=='H') { mass = 1.0;  anum = 1; }
    else if(name[i]=='P') { mass = 31.0; anum = 15; }
    else if(name[i]=='F') { mass = 19.0; anum = 9; }
    else if(name[i]=='B') { mass = 11.0; anum = 5; }
         else { mass = 1.0; anum = 0; } //default
}

}

unique_ptr<Mol_file> Mol_file::recognize(string fname){
    std::string ext = fname.substr(fname.find_last_of(".") + 1);

         if(ext=="xtc")     return unique_ptr<Mol_file>(new XTC_file(fname));
    else if(ext=="trr")     return unique_ptr<Mol_file>(new TRR_file(fname));
    else if(ext=="pdb")     return unique_ptr<Mol_file>(new PDB_file(fname));
    else if(ext=="gro")     return unique_ptr<Mol_file>(new GRO_file(fname));
    else if(ext=="dcd")     return unique_ptr<Mol_file>(new DCD_file(fname));
    else if(ext=="mol2")    return unique_ptr<Mol_file>(new MOL2_file(fname));
    else if(ext=="xyz")     return unique_ptr<Mol_file>(new XYZ_file(fname));
#ifdef USE_TNGIO
    else if(ext=="tng")     return unique_ptr<Mol_file>(new TNG_file(fname));
#endif
#ifdef USE_GROMACS
    else if(ext=="tpr")     return unique_ptr<Mol_file>(new TPR_file(fname));
#endif
#ifdef USE_OPENBABEL
    else if(ext=="pdbqt")   return unique_ptr<Mol_file>(new PDBQT_file(fname));
#endif
    else throw Pteros_error("File extension '{}' is not recognized!",ext);
}

std::unique_ptr<Mol_file> Mol_file::open(string fname, char open_mode)
{
    auto handle = recognize(fname);
    handle->open(open_mode);
    return handle;
}
