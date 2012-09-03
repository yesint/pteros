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

#include "vmd_molfile_plugin_wrapper.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/pdb_cryst.h"
#include <Eigen/Core>

// General molfile_plugin includes
#include "libmolfile_plugin.h"
#include "molfile_plugin.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

/*
    // Internal stuff
    void* file_handle; // Returned by open_file

    // Low-level stuff from VMD molfile_plugin
    std::vector<std::string> extensions;
    void* open_file_read(const char *filepath, const char *filetype, int *natoms);
    int read_structure(void *mydata, int *optflags, molfile_atom_t *atoms);
    int read_bonds(void *v, int *nbonds, int **fromptr, int **toptr,
                          float ** bondorder,int **bondtype,
                          int *nbondtypes, char ***bondtypename);
    int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts);
    void close_read(void *v);

    void *open_file_write(const char *filename, const char *filetype, int natoms);
    int write_structure(void *v, int optflags, const molfile_atom_t *atoms);
    int write_timestep(void *v, const molfile_timestep_t *ts);
    void close_file_write(void *v);
    int read_molecule_metadata(void *v, molfile_metadata_t **metadata);

#include <iostream>
#include <vector>
#include "molfile_plugin.h"
#include "molfile_plugin.h"
#include "readpdb.h"
*/

//extern molfile_plugin_t pdb_plugin;

void box_from_vmd_rep(float fa, float fb, float fc,
                              float alpha, float beta, float gamma, Eigen::Matrix3f& box){
#define XX 0
#define YY 1
#define ZZ 2
#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

    double cosa,cosb,cosg,sing;
    box.fill(0.0);
    box(XX,XX) = fa;

    if ((alpha!=90.0) || (beta!=90.0) || (gamma!=90.0)) {
      if (alpha != 90.0) {
    cosa = cos(alpha*DEG2RAD);
      } else {
    cosa = 0;
      }
      if (beta != 90.0) {
    cosb = cos(beta*DEG2RAD);
      } else {
    cosb = 0;
      }
      if (gamma != 90.0) {
    cosg = cos(gamma*DEG2RAD);
    sing = sin(gamma*DEG2RAD);
      } else {
    cosg = 0;
    sing = 1;
      }
      box(YY,XX) = fb*cosg;
      box(YY,YY) = fb*sing;
      box(ZZ,XX) = fc*cosb;
      box(ZZ,YY) = fc*(cosa - cosb*cosg)/sing;
      box(ZZ,ZZ) = sqrt(fc*fc
             - box(ZZ,XX)*box(ZZ,XX) - box(ZZ,YY)*box(ZZ,YY));
    } else {
      box(YY,YY) = fb;
      box(ZZ,ZZ) = fc;
    }

    // We obtained box as a set of row-vectors. Transform it to column vectors
    box.transposeInPlace();
    box /= 10.0; // Convert to nm
}

VMD_molfile_plugin_wrapper::VMD_molfile_plugin_wrapper(string fname, char open_mode):
    Mol_file(fname,open_mode){    
    handle = NULL;
    w_handle = NULL;
}

VMD_molfile_plugin_wrapper::~VMD_molfile_plugin_wrapper(){
    if(mode=='r'){
        if(handle){
            plugin->close_file_read(handle);
            cout << "Closing file in mode r" << endl;
            handle = NULL;
        }

    } else {
        if(w_handle){
            plugin->close_file_write(w_handle);
            cout << "Closing file in mode w" << endl;
            w_handle = NULL;
        }
    }
}

void VMD_molfile_plugin_wrapper::open(string fname, char open_mode){
    mode = open_mode;
    FILE_FORMATS fmt = recognize_format(fname);    
    if(fmt==accepted_format){

        if(mode=='r'){
            if(handle) throw Pteros_error("Can't open file for reading twice - handle busy!");
            cout << "Opening file '" << fname << "' for reading "
                 << "using VMD plugin '" << plugin->name << "'..." << endl;
            handle = NULL;
            handle = plugin->open_file_read(fname.c_str(), &open_mode, &natoms);
            if(!handle) throw Pteros_error("Can't open file '"+fname + "'!");
            cout << "Number of atoms: " <<natoms << endl;
        } else {
            if(w_handle) throw Pteros_error("Can't open file for writing twice - handle busy!");
            cout << "Opening file '" << fname << "' for writing "
                 << "using VMD plugin '" << plugin->name << "'..." << endl;
            w_handle = NULL;
            // For writing we don't know number of coordinates in selection here
            // so we remember file name and later call open_file_write
            // in write_frame
            stored_write_name = fname;
        }
    } else {
        throw Pteros_error("Format of file "+fname+" is not compatible with selected IO class!");
    }
}

bool VMD_molfile_plugin_wrapper::do_read(System *sys, Frame *frame, Mol_file_content what){

    if(what.structure){
        // READ STRUCTURE:

        if(sys->num_atoms()>0)
            throw Pteros_error("Can't read structure to the system, which is not empty!");

        int flags;
        vector<molfile_atom_t> atoms(natoms);
        plugin->read_structure(handle,&flags,(molfile_atom_t*)&atoms.front());
        // Allocate atoms in the system
        allocate_atoms_in_system(*sys,natoms);
        // Copy atoms to the system
        Atom at;
        for(int i=0; i<natoms; ++i){
            at.name = atoms[i].name;
            at.resname = atoms[i].resname;
            at.resid = atoms[i].resid;
            at.chain = atoms[i].chain[0];
            at.occupancy = atoms[i].occupancy;
            at.beta = atoms[i].bfactor;
            at.mass = atoms[i].mass;

            set_atom_in_system(*sys,i,at);
        }
        sys->assign_resindex();
    }

    if(what.coordinates || what.trajectory){
        // READ FRAME:
        molfile_timestep_t ts;
        //vector<float> buffer(natoms*3);
        //ts.coords = &buffer.front();

        frame->coord.resize(natoms);
        ts.coords = (float*)&frame->coord.front();

        int ret = plugin->read_next_timestep(handle,natoms,&ts);

        if(ret!=MOLFILE_SUCCESS){
            return false;
        }

        for(int i=0; i<natoms; ++i) frame->coord[i].array() /= 10.0;


        // Convert box to our format
        box_from_vmd_rep(ts.A,ts.B,ts.C,ts.alpha,ts.beta,ts.gamma,frame->box);

        frame->t = ts.physical_time;

        return true;

        // Copy coordinates from buffer to our frame and convert to nm
        /*
        frame.coord.resize(natoms);
        int k = 0;
        for(int i=0; i<natoms*3; i+=3){
            frame.coord[k](0) = ts.coords[i]/10.0;
            frame.coord[k](1) = ts.coords[i+1]/10.0;
            frame.coord[k](2) = ts.coords[i+2]/10.0;
            //cout << k << " " <<  frame.coord[k].transpose() << endl;
            ++k;
        }
        */
    }       
}

void VMD_molfile_plugin_wrapper::do_write(Selection &sel, Mol_file_content what){

    if(what.structure){
        // WRITE STRUCTURE:
        if(!w_handle)
            w_handle = plugin->open_file_write(stored_write_name.c_str(), plugin->name, sel.size());

        vector<molfile_atom_t> atoms(sel.size());
        for(int i=0; i<sel.size(); ++i){
            strcpy( atoms[i].name, sel.Name(i).c_str() );
            strcpy( atoms[i].resname, sel.Resname(i).c_str() );
            atoms[i].resid = sel.Resid(i);
            strcpy( atoms[i].chain, &sel.Chain(i) );
            atoms[i].occupancy = sel.Occupancy(i);
            atoms[i].bfactor = sel.Beta(i);
            atoms[i].mass = sel.Mass(i);
        }

        int flags = MOLFILE_OCCUPANCY | MOLFILE_BFACTOR;
        plugin->write_structure(w_handle,flags,&atoms.front());
    }

    if(what.coordinates || what.trajectory){
        // WRITE COORDINATES:
        if(!w_handle)
            w_handle = plugin->open_file_write(stored_write_name.c_str(), plugin->name, sel.size());

        molfile_timestep_t ts;
        int n = sel.size();
        vector<float> buffer(n*3);
        int k = 0;
        for(int i=0; i<n; ++i){
            buffer[k] = sel.X(i)*10.0;
            buffer[k+1] = sel.Y(i)*10.0;
            buffer[k+2] = sel.Z(i)*10.0;
            //cout << k << " " <<  buffer[k] << " " << buffer[k+1] << " " << buffer[k+2] << endl;
            k+=3;
        }
        ts.coords = &buffer.front();
        Eigen::Vector3f v,a;
        box_to_vectors_angles(sel.get_system()->Box(sel.get_frame()),v,a);
        ts.A = v(0)*10.0;
        ts.B = v(1)*10.0;
        ts.C = v(2)*10.0;
        ts.alpha = a(0);
        ts.beta = a(1);
        ts.gamma = a(2);

        ts.physical_time = sel.get_system()->Time(sel.get_frame());

        plugin->write_timestep(w_handle, &ts);
    }
}

