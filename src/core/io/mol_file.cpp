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

#include "pteros/core/mol_file.h"
#include "pteros/core/pteros_error.h"

#include "pdb_file.h"
#include "dcd_file.h"
#include "gro_file.h"
#include "trr_file.h"
#include "xtc_file.h"
#include "pttop_file.h"
#include "tng_file.h"
#include "mol2_file.h"
#include "xyz_file.h"
#include "tpr_file.h"

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


void Mol_file::allocate_atoms_in_system(System &sys, int n){
    sys.atoms.resize(n);
}

void Mol_file::set_atom_in_system(System &sys, int i, Atom &at){
    sys.atoms[i] = at;
}

Atom &Mol_file::atom_in_system(System &sys, int i){
    return sys.atoms[i];
}

void Mol_file::append_atom_in_system(System &sys, Atom &at){
    sys.atoms.push_back(at);
}

Force_field &Mol_file::ff_in_system(System &sys){
    return sys.force_field;
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

float pteros::get_mass_from_atom_name(string& name){
    // Find first character, which is not digit to account for cases like 21C2
    int i = name.find_first_not_of("1234567890");

    if(name[i]=='C') return 12.0;
    else if(name[i]=='O') return 15.0;
    else if(name[i]=='N') return 14.0;
    else if(name[i]=='S') return 32.0;
    else if(name[i]=='H') return 1.0;
    else if(name[i]=='P') return 31.0;
    else return 1.0; //default
}

unique_ptr<Mol_file> Mol_file::recognize(string fname){
    std::string ext = fname.substr(fname.find_last_of(".") + 1);

         if(ext=="xtc")     return unique_ptr<Mol_file>(new XTC_file(fname));
    else if(ext=="trr")     return unique_ptr<Mol_file>(new TRR_file(fname));
    else if(ext=="pdb")     return unique_ptr<Mol_file>(new PDB_file(fname));
    else if(ext=="gro")     return unique_ptr<Mol_file>(new GRO_file(fname));
    else if(ext=="dcd")     return unique_ptr<Mol_file>(new DCD_file(fname));
    else if(ext=="pttop")   return unique_ptr<Mol_file>(new PTTOP_file(fname));
    else if(ext=="tng")     return unique_ptr<Mol_file>(new TNG_file(fname));
    else if(ext=="mol2")    return unique_ptr<Mol_file>(new MOL2_file(fname));
    else if(ext=="xyz")     return unique_ptr<Mol_file>(new XYZ_file(fname));
    else if(ext=="tpr")     return unique_ptr<Mol_file>(new TPR_file(fname));
    else throw Pteros_error("File extension '{}' not recognized!",ext);
}

std::unique_ptr<Mol_file> Mol_file::open(string fname, char open_mode)
{
    auto handle = recognize(fname);
    handle->open(open_mode);
    return handle;
}
