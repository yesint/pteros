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

#include "vmd_molfile_plugin_wrapper.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "pteros/core/utilities.h"
#include "../molfile_plugins/periodic_table.h"
#include <Eigen/Core>
#include <cmath>
// General molfile_plugin includes
#include "molfile_plugin.h"
#include "system_builder.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

using namespace std;
using namespace pteros;
using namespace Eigen;

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



VmdMolfilePluginWrapper::VmdMolfilePluginWrapper(string& fname): FileHandler(fname),
    r_handle(nullptr), w_handle(nullptr)
{

}


void VmdMolfilePluginWrapper::open(char open_mode){
    mode = open_mode;

    if(mode=='r'){
        if(r_handle) throw PterosError("Can't open file for reading twice - handle busy!");
        r_handle = NULL;
        r_handle = plugin->open_file_read(fname.c_str(), &open_mode, &natoms);
        if(!r_handle) throw PterosError("Can't open file '{}'!",fname);
    } else {
        if(w_handle) throw PterosError("Can't open file for writing twice - handle busy!");        
        w_handle = NULL;
    }

}

void VmdMolfilePluginWrapper::close()
{
    if(mode=='r'){
        if(r_handle){
            plugin->close_file_read(r_handle);
            r_handle = NULL;
        }

    } else {
        if(w_handle){
            plugin->close_file_write(w_handle);
            w_handle = NULL;
        }
    }
}

bool VmdMolfilePluginWrapper::do_read(System *sys, Frame *frame, const FileContent &what){

    if(what.atoms()){
        // READ STRUCTURE:
        SystemBuilder builder(sys);

        if(sys->num_atoms()>0)
            throw PterosError("Can't read structure to the system, which is not empty!");

        int flags;
        vector<molfile_atom_t> atoms(natoms);
        plugin->read_structure(r_handle,&flags,(molfile_atom_t*)&atoms.front());
        // Allocate atoms in the system
        builder.allocate_atoms(natoms);
        // Copy atoms to the system
        Atom at;
        int pos,idx;
        for(int i=0; i<natoms; ++i){
            at.name = atoms[i].name;
            at.resname = atoms[i].resname;
            at.resid = atoms[i].resid;
            at.chain = atoms[i].chain[0];
            at.occupancy = atoms[i].occupancy;
            at.beta = atoms[i].bfactor;
            at.charge = atoms[i].charge;
            at.atomic_number = atoms[i].atomicnumber;

            // pdb_plugin guesses mass based on element record, which is
            // often absent. So it is likely that mass will be zero here!
            // Check and guess ourself if there is no mass!
            // We can't use functions from VMD here since atom CA will be
            // treated as calcium for example. Our technique is more primitive
            // and only recognizes few most common elements
            // In this case we also have to guess element itself
            if(atoms[i].mass>0){
                at.mass = atoms[i].mass;
            } else {
                // Guess atomfrom atom name
                get_element_from_atom_name(at.name, at.atomic_number, at.mass);
            }

            builder.set_atom(i,at);
        }

        // Here builder goes out of scope and finalizes the system in destructor
    }

    if(what.atoms()) sys->assign_resindex();

    if(what.coord() || what.traj()){
        // READ FRAME:
        molfile_timestep_t ts;        
        // Set zeros to box variables
        ts.A = ts.B = ts.C = ts.alpha = ts.beta = ts.gamma = 0.0;

        frame->coord.resize(natoms);
        ts.coords = (float*)&frame->coord.front();

        int ret = plugin->read_next_timestep(r_handle,natoms,&ts);

        if(ret!=MOLFILE_SUCCESS){
            return false;
        }

        for(int i=0; i<natoms; ++i) frame->coord[i].array() /= 10.0;

        // Convert box to our format
        Matrix3f b;
        b.fill(0.0);
        if(ts.A*ts.B*ts.C){
            // Only convert if all three vectors are non-zero
            box_from_vmd_rep(ts.A,ts.B,ts.C,ts.alpha,ts.beta,ts.gamma,b);
        }
        frame->box.set_matrix(b);

        frame->time = ts.physical_time;

        return true;
    }

    // If we are here than something is wrong
    return false;
}

void VmdMolfilePluginWrapper::do_write(const Selection &sel, const FileContent &what) {

    if(what.atoms()){
        // WRITE STRUCTURE:        
        if(!w_handle)
            w_handle = plugin->open_file_write(fname.c_str(), plugin->name, sel.size());

        vector<molfile_atom_t> atoms(sel.size());        
        for(int i=0; i<sel.size(); ++i){
            strcpy( atoms[i].name, sel.name(i).c_str() );
            strcpy( atoms[i].resname, sel.resname(i).c_str() );
            atoms[i].resid = sel.resid(i);
            //stringstream ss;
            //ss << sel.chain(i);
            //strcpy( atoms[i].chain, ss.str().c_str() );
            atoms[i].chain[0] = sel.chain(i);
            atoms[i].chain[1] = '\0';
            atoms[i].occupancy = sel.occupancy(i);
            atoms[i].bfactor = sel.beta(i);
            atoms[i].mass = sel.mass(i);
            atoms[i].charge = sel.charge(i);
            atoms[i].atomicnumber = sel.atomic_number(i);

            // For MOL2 we also need to set atom type as a string
            strcpy( atoms[i].type, sel.element_name(i).c_str() );
        }
        int flags = MOLFILE_OCCUPANCY | MOLFILE_BFACTOR | MOLFILE_ATOMICNUMBER
                    | MOLFILE_CHARGE | MOLFILE_MASS;
        plugin->write_structure(w_handle,flags,&atoms.front());        
    }

    if(what.coord() || what.traj()){
        // WRITE COORDINATES:
        if(!w_handle)
            w_handle = plugin->open_file_write(fname.c_str(), plugin->name, sel.size());

        molfile_timestep_t ts;
        int n = sel.size();
        vector<float> buffer(n*3);
        int k = 0;
        for(int i=0; i<n; ++i){
            buffer[k] = sel.x(i)*10.0;
            buffer[k+1] = sel.y(i)*10.0;
            buffer[k+2] = sel.z(i)*10.0;
            k+=3;
        }
        ts.coords = &buffer.front();
        ts.velocities = NULL; // No velocities currently supported
        // Only convert periodic box if it is present
        if(sel.box().is_periodic()){
            Eigen::Vector3f v,a;
            sel.box().to_vectors_angles(v,a);
            ts.A = v(0)*10.0;
            ts.B = v(1)*10.0;
            ts.C = v(2)*10.0;
            ts.alpha = a(0);
            ts.beta = a(1);
            ts.gamma = a(2);
        } else {
            ts.A = 0.0;
            ts.B = 0.0;
            ts.C = 0.0;
            ts.alpha = 0;
            ts.beta = 0;
            ts.gamma = 0;
        }

        ts.physical_time = sel.get_system()->time(sel.get_frame());

        plugin->write_timestep(w_handle, &ts);
    }
}

//--------------------------------------------------------
// molfile plugins registration bootstrap
//--------------------------------------------------------

#define IMPORT_PLUGIN(name) \
    VMDPLUGIN_EXTERN int name##plugin_init(); \
    VMDPLUGIN_EXTERN int name##plugin_register(void *v, vmdplugin_register_cb cb); \
    VMDPLUGIN_EXTERN int name##plugin_fini();

#define REGISTER_PLUGIN(name,ret) \
    name##plugin_init(); \
    name##plugin_register(nullptr, &register_cb); \
    ret[cur_name] = cur_plugin;

IMPORT_PLUGIN(pdb)
IMPORT_PLUGIN(dcd)
IMPORT_PLUGIN(xyz)

#ifdef USE_TNGIO
IMPORT_PLUGIN(tng)
#endif

molfile_plugin_t *cur_plugin;
string cur_name;

static int register_cb(void *v, vmdplugin_t *p) {
  cur_name = string(p->name);
  cur_plugin = (molfile_plugin_t *)p;
  return VMDPLUGIN_SUCCESS;
}


std::map<string,molfile_plugin_t*> register_all_plugins(){    
    std::map<string,molfile_plugin_t*> ret;

    REGISTER_PLUGIN(pdb,ret)
    REGISTER_PLUGIN(dcd,ret)
    REGISTER_PLUGIN(xyz,ret)    

#ifdef USE_TNGIO
    REGISTER_PLUGIN(tng,ret)
#endif

    // Debug output on loaded plugins
    LOG()->debug("Registered VMD molfile plugins:");
    for(auto& item: ret){
         LOG()->debug("{}", item.first);
    }

    return ret;
}

std::map<string,molfile_plugin_t*> VmdMolfilePluginWrapper::molfile_plugins = register_all_plugins();




