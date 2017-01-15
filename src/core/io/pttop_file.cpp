/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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

#include "pttop_file.h"
#include "pteros/core/pteros_error.h"
#include <fstream>

using namespace std;
using namespace pteros;
using namespace Eigen;


void PTTOP_file::open(char open_mode)
{
    if(open_mode=='w')
        throw Pteros_error("PTTOP files could be produced by the dedicated script only!");
}

PTTOP_file::~PTTOP_file(){
}

bool PTTOP_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){
    string line;
    stringstream ss;
    int dum;
    float val;

    ifstream f(fname.c_str());
    if(!f) throw Pteros_error("Can't open PTTOP file '"+fname+"'' for reading");

    cout << "Reading Pteros topology file '"<< fname << "'..." << endl;

    // coulomb_type
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).coulomb_type;

    // coulomb_modifier
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).coulomb_modifier;

    // rcoulomb
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).rcoulomb;

    // epsilon_r
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).epsilon_r;

    // epsilon_rf
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).epsilon_rf;

    // rcoulomb_switch
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).rcoulomb_switch;

    // vdw_type
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).vdw_type;

    // vdw_modifier
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).vdw_modifier;

    // rvdw_switch
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).rvdw_switch;

    // rvdw
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    if(what.top()) ss >> ff_in_system(*sys).rvdw;

    // Number of atoms
    int natoms;
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    ss >> natoms;
    if(frame) frame->coord.resize(natoms);
    getline(f,line); // Comment
    // Read atoms
    for(int i=0;i<natoms;++i){        
        Atom at;
        getline(f,line);
        ss.clear(); ss.str(line);
        ss >> dum  >> at.name >> at.type_name >> at.type >> at.resid
           >> at.resindex >> at.resname >> at.mass >> at.charge;
        if(what.coord()){
           ss >> frame->coord[i](0) >> frame->coord[i](1) >> frame->coord[i](2);
        }
        if(what.atoms()){
            // Add atoms to system
            append_atom_in_system(*sys,at);
        } else {
            // Update atoms in system
            atom_in_system(*sys,i) = at;
        }
    }
    cout << "\tRead " << natoms << " atoms" << endl;

    // Bonds
    int nbonds;
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    ss >> nbonds;    
    ff_in_system(*sys).bonds.resize(nbonds);    
    getline(f,line); // Comment
    for(int i=0;i<nbonds;++i){        
        getline(f,line);
        ss.clear(); ss.str(line);
        ss >> ff_in_system(*sys).bonds[i](0) >> ff_in_system(*sys).bonds[i](1);
    }
    cout << "\tRead " << nbonds << " bonds" << endl;

    // Box
    getline(f,line); // Comment
    Matrix3f b;
    for(int i=0; i<3; ++i){
        getline(f,line);
        ss.clear(); ss.str(line);
        if(what.coord()){
            ss >> b(i,0) >> b(i,1) >> b(i,2);
        }
    }
    if(frame) frame->box.modify(b);

    // All the rest should be read only if topology is requested
    if(!what.top()) return true;

    // Number of charge groups
    int ncg;
    getline(f,line); // Coment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    ss >> ncg;
    ff_in_system(*sys).charge_groups.resize(ncg);
    getline(f,line); // Coment
    // Read charge groups
    for(int i=0;i<ncg;++i){
        getline(f,line);
        ss.clear(); ss.str(line);
        ss >> ff_in_system(*sys).charge_groups[i](0)
           >> ff_in_system(*sys).charge_groups[i](1);
    }
    cout << "\tRead " << ncg << " charge groups" << endl;

    // number of atoms with exclusions
    int nexcl;
    getline(f,line); // Coment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    ss >> nexcl;
    ff_in_system(*sys).exclusions.resize(natoms);
    getline(f,line); // Coment
    // Read charge groups
    for(int i=0;i<nexcl;++i){
        getline(f,line);
        ss.clear(); ss.str(line);
        ss >> dum; // atom index itself
        while(ss){
            ss >> dum;
            ff_in_system(*sys).exclusions[i].insert(dum);
        }
    }
    cout << "\tRead " << nexcl << " exclusion groups" << endl;

    // Size of LJ matrix
    int nLJ;
    getline(f,line); // Comment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    ss >> nLJ;
    ff_in_system(*sys).LJ_C6.resize(nLJ,nLJ);
    ff_in_system(*sys).LJ_C12.resize(nLJ,nLJ);
    getline(f,line); // Comment
    // Read C6 matrix
    for(int i=0;i<nLJ;++i){
        getline(f,line);
        ss.clear(); ss.str(line);
        for(int j=0;j<nLJ;++j){
            ss >> val;
            ff_in_system(*sys).LJ_C6(i,j) = val;
        }
    }
    getline(f,line); // Comment
    // Read C12 matrix
    for(int i=0;i<nLJ;++i){
        getline(f,line);
        ss.clear(); ss.str(line);
        for(int j=0;j<nLJ;++j){
            ss >> val;
            ff_in_system(*sys).LJ_C12(i,j) = val;
        }
    }
    cout << "\tRead LJ interaction matrices of size " << nLJ << endl;

    // fudgeQQ
    getline(f,line); // Coment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    ss >> ff_in_system(*sys).fudgeQQ;

    // number of LJ14 interaction types
    int nLJ14types;
    getline(f,line); // Coment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    ss >> nLJ14types;
    ff_in_system(*sys).LJ14_interactions.resize(nLJ14types);
    getline(f,line); // Coment
    // Read LJ14 types
    for(int i=0;i<nLJ14types;++i){
        getline(f,line);
        ss.clear(); ss.str(line);
        ss >> dum >> ff_in_system(*sys).LJ14_interactions[i](0)
                  >> ff_in_system(*sys).LJ14_interactions[i](1);
    }
    cout << "\tRead " << nLJ14types << " LJ-14 interaction types" << endl;

    // number of LJ14 pairs
    int nLJ14pairs;
    getline(f,line); // Coment
    getline(f,line); // Value
    ss.clear(); ss.str(line);
    ss >> nLJ14pairs;
    getline(f,line); // Coment
    // Read LJ14 types
    for(int i=0;i<nLJ14pairs;++i){
        getline(f,line);
        ss.clear(); ss.str(line);
        int a,b,n;
        ss >> a >> b >> n;
        ff_in_system(*sys).LJ14_pairs.insert( make_pair(a*nLJ14types+b,n) );
    }
    cout << "\tRead " << nLJ14pairs << " LJ-14 atom pairs" << endl;

    ff_in_system(*sys).setup_kernels(); // Setup kernels
    ff_in_system(*sys).ready = true; // ff is ready to use

    return true;
}
