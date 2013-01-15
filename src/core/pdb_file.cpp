/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
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

#include "pdb_file.h"
#include "molfile_plugin.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/pdb_cryst.h"
#include <iomanip>

using namespace std;
using namespace pteros;
using namespace Eigen;

#ifdef VMD_PDB

extern molfile_plugin_t pdb_plugin;

PDB_file::PDB_file(string fname, char open_mode): VMD_molfile_plugin_wrapper(fname, open_mode)
{
    plugin = &pdb_plugin;
    accepted_format = PDB_FILE;
    open(fname,open_mode);
}


#else

PDB_file::PDB_file(string fname, char open_mode): Mol_file(fname, open_mode)
{
    if(open_mode=='r'){
        f.open(fname.c_str(),ios_base::in);
        if(!f) throw Pteros_error("Can't open PDB file '"+fname+"'' for reading");
    } else {
        f.open(fname.c_str(),ios_base::out);
        if(!f) throw Pteros_error("Can't open PDB file '"+fname+"'' for writing");
    }
}

PDB_file::~PDB_file(){
    if(f) f.close();
}


bool PDB_file::do_read(System *sys, Frame *frame, Mol_file_content what){

    string line;
    // Temporary Atom object
    Atom tmp_atom;
    // Temporary coordinates
    Vector3f tmp_coor;
    // Temporary index
    long tmp_index;
    long index = 0;
    stringstream ss;
    string dum;
    char ch;

    float box_x, box_y, box_z, box_a, box_b, box_c;

    frame->coord.clear();

    //Start reading file
    while(getline(f,line)){
        ss.clear();
        ss.str(line);
        ss >> dum;
        if(dum=="ATOM" || dum=="HETATM"){
            ss >> dum; //index. we don't need it
            ss >> tmp_atom.name;
            ss >> tmp_atom.resname;
            // Get chain at character level (may be space)
            ss >> noskipws >> ch >> tmp_atom.chain >> skipws;
            //  If the chain is space, then change it to X
            if(tmp_atom.chain==' ') tmp_atom.chain = 'X';
            ss >> tmp_atom.resid;
            ss >> tmp_coor[0] >> tmp_coor[1] >> tmp_coor[2];
            ss >> tmp_atom.occupancy;
            ss >> tmp_atom.beta;
            ss >> setw(4) >> tmp_atom.tag;

            if(what.structure){ // Reading structure
                // Assign masses
                tmp_atom.mass = get_mass_from_atom_name(tmp_atom.name);
                tmp_atom.type = -1; //Undefined type so far
                // Add new atom to the system
                append_atom_in_system(*sys,tmp_atom);
            }

            if(what.coordinates){
                tmp_coor /= 10.0; // Convert to nm
                // Add column of coordinates
                frame->coord.push_back(tmp_coor);
            }

            // Go to next atom
            index++;
        } else if(dum=="CRYST1" && what.structure) {
            // Converting the box entry is rather tedious.
            // Adapted routine from GROMACS 4.x is used to avoid errors
            read_pdb_box(line.c_str(), frame->box);
            frame->box /= 10.0; //Convert to nm
        }
    } //Done reading

    return true;
}


void PDB_file::do_write(Selection &sel, Mol_file_content what){
    int n = sel.size();    

    if(what.structure && what.coordinates){
        // Write periodic box
        f << write_pdb_box(sel.get_system()->Box(sel.get_frame())*10.0); //Convert to Angstroms

        for(int i=0;i<n;i++){
            f << "ATOM  " << right
              << setw(5) << i+1
              << setw(4) << sel.Name(i)
              << setw(5) << sel.Resname(i)
              << setw(2) << sel.Chain(i)
              << setw(4) << sel.Resid(i)
              << fixed << setprecision(3) << setw(12) << sel.X(i)*10.0 //Convert back to Angstorms
              << setw(8) << sel.Y(i)*10.0 << setw(8) << sel.Z(i)*10.0
              << setprecision(2) << setw(6) << sel.Occupancy(i)
              << setw(6) << sel.Beta(i)
              << setw(10) << sel.Tag(i) << endl;
        }
    } else {
        throw Pteros_error("Atoms and coordinates should always be written to PDB file!");
    }
}


#endif
