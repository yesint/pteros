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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/grid_search.h"
#include "pteros/core/pdb_cryst.h"
#include "pteros/core/format_recognition.h"

#include "pteros/core/mol_file.h"
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace pteros;
using namespace Eigen;

// Base constructor of the system class
System::System() {

}

// Construnt system from file
System::System(string fname) {
    //read_structure(fname);
    clear();
    load(fname);
}

System::System(const System& other){
    clear(true);
    atoms = other.atoms;
    traj = other.traj;
}

System& System::operator=(System other){
    clear(true);
    atoms = other.atoms;
    traj = other.traj;
    return *this;
}

// Clear the system (i.e. before reading new system from file)
void System::clear(bool delete_selections){
    atoms.clear();
    traj.clear();
    if(delete_selections) notify_signal(SYSTEM_CLEARED,0,0);
}

void System::check_num_atoms_in_last_frame(){
    if(Frame_data(num_frames()-1).coord.size()!=num_atoms())
        throw Pteros_error("File contains "
                           +boost::lexical_cast<string>(Frame_data(num_frames()-1).coord.size())
                           +" atoms while the system has "
                           +boost::lexical_cast<string>(num_atoms())
                           );
}

// Load structure or trajectory
void System::load(string fname, int b, int e, int skip){
    // Create an IO file reader
    boost::shared_ptr<Mol_file> f = io_factory(fname,'r');

    int num_stored = 0;    
    // Do we have some structure?
    if(num_atoms()>0){
        // We have atoms already, so read only coordinates
        if(!f->get_content_type().coordinates && !f->get_content_type().trajectory)
            throw Pteros_error("File reader for file '"+fname
                               +"' is not capable of appending frames to the system!");

        // Check if we can read multiple coordinates
        if(f->get_content_type().trajectory){

            Mol_file_content c;
            c.trajectory = true;

            // Sanity check for frame range
            if((e<b && e!=-1)|| b<0)
                throw Pteros_error("Invalid frame range for reading!");

            int cur = 0; // This holds real frame index in trajectory

            // Skip frames if needed
            if(b>0){
                cout << "Skipping " << b << " frames..." << endl;
                Frame skip_fr;
                for(int i=0;i<b;++i){
                    f->read(NULL,&skip_fr,c);
                    cur++;
                }
            }

            int first = num_frames(); // Remember start

            cout << "Reading..."<<endl;

            int actually_read = 0;            

            while(true){
                // End frame reached?
                if(cur==e && e!=-1) break;

                // Append new frame where the data will go
                Frame fr;
                frame_append(fr);
                // Try to read into it
                bool ok = f->read(NULL,&Frame_data(num_frames()-1),c);
                if(!ok){
                    frame_delete(num_frames()-1); // Remove last frame - it's invalid
                    break;
                }

                check_num_atoms_in_last_frame();

                ++cur;
                ++actually_read;

                if(skip>0 && actually_read%skip!=0){
                    frame_delete(num_frames()-1); // Remove last frame - it's invalid
                    continue;
                } else {
                    actually_read = 0;
                }

                // If we are here new frame is already accepted
                ++num_stored;
            }        
        } else if(f->get_content_type().coordinates) {
            Mol_file_content c;
            c.coordinates = true;
            // File contains single frame
            // Append new frame where the data will go
            Frame fr;
            frame_append(fr);
            // Read it
            f->read(NULL,&Frame_data(num_frames()-1),c);
            check_num_atoms_in_last_frame();
            ++num_stored;
        }
    } else {
        // We don't have atoms yet, so we will read everything possible
        // Append new frame where the data will go
        Frame fr;
        frame_append(fr);
        Mol_file_content c = f->get_content_type();
        f->read(this,&Frame_data(num_frames()-1),c);

        check_num_atoms_in_last_frame();
        ++num_stored;

        assign_resindex();
    }

    cout << "Stored " << num_stored << " frames. Now " << num_frames() << " frames in the System" << endl;
}

/*
void System::compute_bonds(){
    Selection s(*this,"all");
    Grid_searcher(0.16,s,bonds);
}
*/

// Update all selections
void System::update_selections(){
    notify_signal(TOPOLOGY_CHANGED,0,0);
}

// Destructor of the system class
System::~System() {}

void System::set_frame(int fr){
    if(fr<0 || fr>traj.size()) throw Pteros_error("Invalid frame!");
    notify_signal(FRAME_CHANGE_REQUESTED,fr,0);
}

void System::frame_dup(int fr){
    if(fr<0 || fr>=traj.size())
    	throw Pteros_error("Invalid frame for duplication!");
    traj.push_back(traj[fr]);
}

void System::frame_copy(int fr1, int fr2){
    if(fr1<0 || fr1>=traj.size() || fr2<0 || fr2>=traj.size())
    	throw Pteros_error("Invalid frame for copying!");
    traj[fr2] = traj[fr1];
    // Coordinate-dependent selections, which point to fr2 should be updated
    notify_signal(COORDINATES_CHANGED,fr2,fr2);
}

// Delete the range of frames. e = -1 is default
void System::frame_delete(int b, int e){
    int i;    

    if(e==-1) e = num_frames()-1;
    if(e<b || b<0 || e>num_frames()-1) throw Pteros_error("Invalid frame range for deletion");    

    // Get iterators for deletion
    vector<Frame>::iterator b_it, e_it;
    b_it = traj.begin();
    for(i=0;i<b;++i) b_it++;
    e_it = b_it;
    for(;i<e;++i) e_it++;
    e_it++; //Go one past the end

    traj.erase(b_it,e_it);

    // Check if there are some frames left. If not print the warning
    // that all selections are invalid!
    if(traj.size()==0) cout << "All frames are deleted. All selections are now INVALID!";

    // Selections, which point to frames >b now become invalid
    // Make them point to frame zero
    notify_signal(FRAMES_DELETED,b,e);
}

void System::get_box_vectors_angles(int fr, Vector3f& vectors, Vector3f& angles) const {
    // Call Gromacs function from pdb_cryst
    box_to_vectors_angles(traj[fr].box, vectors, angles);
}

bool System::is_box_triclinic() const{
    return (traj[0].box(0,1) || traj[0].box(0,2)
            || traj[0].box(1,0) || traj[0].box(1,2)
            || traj[0].box(2,0) || traj[0].box(2,1) );
}

void System::frame_append(const Frame& fr){
    traj.push_back(fr);
}

Frame& System::Frame_data(int fr) {
    //if(fr<0 || fr>=traj.size())
    //	throw Pteros_error("Invalid frame!");
    return traj[fr];
}

void System::assign_resindex(){
    //cout << "Assigning resindex..." << endl;
    int curres = atoms[0].resid;
    int curchain = atoms[0].chain;
    int cur = 0;
    for(int i=0; i<atoms.size(); ++i){
        if( atoms[i].resid!=curres || atoms[i].chain!=curchain ){
            ++cur;
            curres = atoms[i].resid;
            curchain = atoms[i].chain;
        }
        atoms[i].resindex = cur;
    }
}

void System::atoms_dup(const vector<int>& ind, Selection* res_sel){
    // Sanity check
    if(!ind.size()) throw Pteros_error("No atoms to duplicate!");
    for(int i=0; i<ind.size(); ++i){
        if(ind[i]<0 || ind[i]>atoms.size()-1)
            throw Pteros_error("Invalid index for atom duplication!");
    }

    // Duplicate atoms
    int first_added = atoms.size();
    int last_added = atoms.size()+ind.size()-1;
    // Prepare by increasing capacity of vectors
    atoms.reserve(atoms.size()+ind.size());
    for(int j=0; j<traj.size(); ++j){
        traj[j].coord.reserve(atoms.size()+ind.size());
    }

    // Now add atoms
    for(int i=0; i<ind.size(); ++i){
        // Add new atom
        atoms.push_back(atoms[ind[i]]);
        // Add new coordinate slot
        for(int j=0; j<traj.size(); ++j){
            traj[j].coord.push_back(traj[j].coord[ind[i]]);
        }
    }

    if(res_sel) res_sel->modify(*this,first_added,last_added);
}

void System::atoms_add(const vector<Atom>& atm, const vector<Vector3f>& crd, Selection* res_sel){
    // Sanity check
    if(!atm.size()) throw Pteros_error("No atoms to add!");
    if(atm.size()!=crd.size()) throw Pteros_error("Wrong number of coordinates for adding atoms!");

    int first_added = atoms.size();
    int last_added = atoms.size()+atm.size()-1;
    // Prepare by increasing capacity of vectors
    atoms.reserve(atoms.size()+atm.size());
    for(int j=0; j<traj.size(); ++j){
        traj[j].coord.reserve(atoms.size()+atm.size());
    }
    // Now add atoms
    for(int i=0; i<atm.size(); ++i){
        // Add new atom
        atoms.push_back(atm[i]);
        // Add new coordinate slot
        for(int j=0; j<traj.size(); ++j){
            traj[j].coord.push_back(crd[i]);
        }
    }

    if(res_sel) res_sel->modify(*this,first_added,last_added);
}

void System::atoms_delete(const std::vector<int> &ind){
    int i,fr;

    // Sanity check
    if(!ind.size()) throw Pteros_error("No atoms to delete!");
    for(int i=0; i<ind.size(); ++i){
        if(ind[i]<0 || ind[i]>atoms.size()-1)
            throw Pteros_error("Invalid index for atom deletion!");
    }

    // Mark atoms for deletion by assigning negative mass
    for(i=0;i<ind.size();++i)
        atoms[ind[i]].mass = -1.0;

    // Cycle over all atoms and keep only those with positive mass
    vector<pteros::Atom> tmp = atoms;
    atoms.clear();
    for(i=0;i<tmp.size();++i){
        if(tmp[i].mass>=0) atoms.push_back(tmp[i]);
    }

    // Now cycle over trajectory and keep only corresponding coordinates
    vector<Vector3f> tmp_coord;
    for(fr=0; fr<num_frames(); ++fr){
        // Make a copy of traj coords
        tmp_coord = traj[fr].coord;
        traj[fr].coord.clear();
        for(i=0;i<tmp.size();++i){
            if(tmp[i].mass>=0) traj[fr].coord.push_back(tmp_coord[i]);
        }
    }
    // Reassign residue indexes
    //system->assign_resindex();
    //system->update_selections();
}

void System::append(const System &sys){
    //Sanity check
    if(num_frames()!=sys.num_frames()) throw Pteros_error("Can't merge systems with different number of frames!");
    // Merge atoms
    copy(sys.atoms.begin(),sys.atoms.end(),back_inserter(atoms));
    // Merge coordinates
    for(int fr=0; fr<num_frames(); ++fr)
        copy(sys.traj[fr].coord.begin(),sys.traj[fr].coord.end(),back_inserter(traj[fr].coord));
    // Reassign resindex
    assign_resindex();
    // Update selections
    update_selections();
}

inline void wrap_coord(Vector3f& point, const Vector3f& box_dim){
    int i;
    float intp,fracp;
    for(i=0;i<3;++i){
        fracp = std::modf(point(i)/box_dim(i),&intp);
        if(fracp<0) fracp = fracp+1;
        point(i) = box_dim(i)*fracp;
    }
}

float System::distance(const Eigen::Vector3f &p1, const Eigen::Vector3f &p2, int fr, bool is_periodic) const {
    if(is_periodic && traj[fr].box.array().abs().sum()>0){ // If box is not set, it will be zero
        // Get box dimension
        Vector3f box_dim = traj[fr].box.colwise().norm();
        Vector3f pp1 = p1, pp2 = p2;
        wrap_coord(pp1,box_dim);
        wrap_coord(pp2,box_dim);
        // For each dimension measure periodic distance
        Vector3f v = (pp2-pp1).array().abs();
        for(int i=0;i<3;++i)
            if(v(i)>0.5*box_dim(i)) v(i) = box_dim(i)-v(i);
        return v.norm();
    } else {
        return (p2-p1).norm();
    }
}

float System::distance(int i, int j, int fr, bool is_periodic) const {
    return distance(traj[fr].coord[i], traj[fr].coord[j], fr, is_periodic);
}

void System::wrap_to_box(int frame, Eigen::Vector3f& point) const {
    int i;
    float intp,fracp;
    // Get box vectors norms
    Vector3f box_dim = traj[frame].box.colwise().norm();
    wrap_coord(point,box_dim);
}

Vector3f System::get_closest_image(Eigen::Vector3f &point, Eigen::Vector3f &target, int fr, bool do_wrapping) const {
    if(traj[fr].box.array().abs().sum()>0){ // If box is not set, it will be zero
        // Get box dimension
        Vector3f box_dim = traj[fr].box.colwise().norm();
        // Wrap point and target
        Vector3f p = point, t = target;
        if(do_wrapping){
            wrap_coord(p,box_dim);
            wrap_coord(t,box_dim);
        }

        Vector3f v = (p-t).array();
        //cout << v(0) <<  " -- " << 0.5*box_dim(0) << endl;
        for(int i=0;i<3;++i)            
            if(abs(v(i))>0.5*box_dim(i)){                
                // Need to translate this dimension
                v(i)>0 ? p(i)-=box_dim(i) : p(i)+=box_dim(i);
            }

        return p;
    } else {
        return point;
    }
}

void System::wrap_all(int fr){
    for(int i=0;i<num_atoms();++i){
        wrap_to_box(fr,XYZ(i,fr));
    }
}

inline float LJ_en_kernel(float sig, float eps, float r){
    float tmp = sig/r;
    tmp = tmp*tmp; // This gets (s/r)^2
    tmp = tmp*tmp*tmp; // This gets (s/r)^6
    return 4.0*eps*(tmp*tmp-tmp);
}

#define ONE_4PI_EPS0      138.935456

inline float Coulomb_en_kernel(float q1, float q2, float r){
    return ONE_4PI_EPS0*q1*q2/r;
}

string Energy_components::to_str(){
    return    boost::lexical_cast<string>(total) + " "
            + boost::lexical_cast<string>(lj_sr) + " "
            + boost::lexical_cast<string>(lj_14) + " "
            + boost::lexical_cast<string>(q_sr) + " "
            + boost::lexical_cast<string>(q_14);
}

void System::add_non_bond_energy(Energy_components &e, int a1, int a2, int frame, bool is_periodic)
{
    // First check if this pair is not in exclusions
    if( force_field.exclusions[a1].count(a2) == 0 ){
        // Required at1 < at2
        int at1,at2;
        if(a1<a2){
            at1 = a1;
            at2 = a2;
        } else {
            at1 = a2;
            at2 = a1;
        }

        int N = force_field.LJ14_interactions.size();
        float r = distance(XYZ(at1,frame),XYZ(at2,frame),frame,is_periodic);
        // Check if this is 1-4 pair
        auto it = force_field.LJ14_pairs.find(at1*N+at2);
        if( it == force_field.LJ14_pairs.end() ){
            // Normal, not 1-4
            e.lj_sr += LJ_en_kernel(force_field.LJ_C6(atoms[at1].type,atoms[at2].type),
                                   force_field.LJ_C12(atoms[at1].type,atoms[at2].type),
                                   r);
            e.q_sr += Coulomb_en_kernel(atoms[at1].charge,
                                       atoms[at2].charge,
                                       r);
            e.total += e.lj_sr + e.q_sr;
        } else {
            // 1-4
            e.lj_14 += LJ_en_kernel(force_field.LJ14_interactions[it->second](0),
                                   force_field.LJ14_interactions[it->second](1),
                                   r);
            e.q_14 += Coulomb_en_kernel(atoms[at1].charge,
                                       atoms[at2].charge,
                                       r)
                    * force_field.fudgeQQ;
            e.total += e.lj_14 + e.q_14;
        }
    }
}
