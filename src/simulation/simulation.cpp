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

#include "pteros/simulation/simulation.h"
#include "pteros/simulation/gromacs_ff.h"
#include "pteros/core/format_recognition.h"
#include "pteros/core/pteros_error.h"
#include <boost/lexical_cast.hpp>

#define ONE_4PI_EPS0      138.935456

using namespace std;
using namespace pteros;
using namespace Eigen;

string Energy_components::to_str(){
    return    boost::lexical_cast<string>(total) + " "
            + boost::lexical_cast<string>(lj_sr) + " "
            + boost::lexical_cast<string>(lj_14) + " "
            + boost::lexical_cast<string>(q_sr) + " "
            + boost::lexical_cast<string>(q_14);
}


Simulation::Simulation(){
    //setup_ok = false; // Not complete yet
}

// Copy constructor
Simulation::Simulation(const Simulation& other){
    system = other.system;
    frame = other.frame;
    Natoms = other.Natoms;
    velocity = other.velocity;
    force = other.force;
    inv_mass = other.inv_mass;
    sigma = other.sigma;
    epsilon = other.epsilon;
    charge = other.charge;
    comb_rule_function = other.comb_rule_function;
    scale14_coulomb = other.scale14_coulomb;
    scale14_lj = other.scale14_lj;
    constr = other.constr;
    angles = other.angles;
    dihedrals = other.dihedrals;
    bond_dist = other.bond_dist;
}

// Assignment operator
Simulation& Simulation::operator=(Simulation other){
    system = other.system;
    frame = other.frame;
    Natoms = other.Natoms;
    velocity = other.velocity;
    force = other.force;
    inv_mass = other.inv_mass;
    sigma = other.sigma;
    epsilon = other.epsilon;
    charge = other.charge;
    comb_rule_function = other.comb_rule_function;
    scale14_coulomb = other.scale14_coulomb;
    scale14_lj = other.scale14_lj;
    constr = other.constr;
    angles = other.angles;
    dihedrals = other.dihedrals;
    bond_dist = other.bond_dist;
    return *this;
}

Simulation::Simulation(System& sys, std::string top_file){
    //setup_ok = false; // Not complete yet
    create(sys,top_file);
    frame = 0;
}

Simulation::~Simulation(){
}

void Simulation::create(System& sys, std::string top_file){
    system = &sys;
    // Set number of atoms
    Natoms = sys.num_atoms();
    // Init fast arrays
    sigma.resize(Natoms);
    epsilon.resize(Natoms);
    charge.resize(Natoms);
    inv_mass.resize(Natoms);
    force.resize(3,Natoms);
    velocity.resize(3,Natoms);

    // Load topology.
    // At some point we'll probably support also charmm PSF format

    FILE_FORMATS fmt = recognize_format(top_file);
    switch(fmt){
        case TOP_FILE:
            read_gromacs_top(top_file);
            break;
        default:
            throw Pteros_error("File "+top_file+" is not valid topology file!");
    }
}

/*
void Simulation::add_excluded_groups(vector<string>& gr){
    if(gr.size()){
        cout << "Preparing excluded groups..." << endl;
        excluded_groups.clear(); // Discard all existing excluded groups
        string non_excl("not (");
        for(int i=0;i<gr.size();++i){
            excluded_groups.push_back(Selection(*all.get_system(),gr[i]));
            cout << "\tAdded excluded group " << excluded_groups.size()
                 <<" of size\t" << excluded_groups.back().size() << endl;

            non_excl += "(" + gr[i] + ") ";
            if(i!=gr.size()-1) non_excl += "or ";
        }
        non_excl += ")";
        // Make selection containing atoms not included in any excluded group
        not_excluded.modify(*all.get_system(),non_excl);
    } else
        not_excluded.modify(*all.get_system(),"all");

    cout << "\tCreated not excluded group of size\t" << not_excluded.size() << endl << endl;
}

void Simulation::setup(){
    // Test validity of setup here
    //...

    // If excluded groups were not provided, create not_excluded here
    if(excluded_groups.size()==0 && not_excluded.size()==0)
        not_excluded.modify(*all.get_system(),"all");

    // Init grid searcher
    //Absolute indexes will be reported
    searcher.create(params.nb_cutoff, *all.get_system(), true, params.is_periodic);

    // Initial velocity
    velocity.fill(0.0);
    // Initial force
    compute_force();

    setup_ok = true;
    cout << "Setup complete:" << endl;
    // Print some statistics
    cout << "\tNumber of atoms:\t" << Natoms << endl;
    cout << "\tNumber of bonds:\t" << constr.size() << endl;
    cout << "\tNumber of angles:\t" << angles.size() << endl;
    cout << "\tNumber of dihedrals:\t" << dihedrals.size() << endl;
    cout << "\tNon-bond cut-off:\t" << params.nb_cutoff << endl;
    cout << "\tLJ combination rule:\t";
    switch(comb_rule){
        case COMB_RULE_GROMOS:
            cout << "GROMOS (GROMACS type 1)";
            break;
        case COMB_RULE_AMBER:
            cout << "AMBER (GROMACS type 2)";
            break;
        case COMB_RULE_OPLS:
            cout << "OPLS (GROMACS type 3)";
            break;
    }
    cout << endl;
    cout << "Ready to run." << endl;
}

void Simulation::update_nlist(){
    // Clear old nlist
    nlist.clear();

    int Ngr = excluded_groups.size();

    // Compute neighbours between all excluded groups
    // and between each excluded and non-excluded
    if(Ngr){
        int i,j;

        for(i=0;i<Ngr-1;++i){
            searcher.search(excluded_groups[i], not_excluded, nlist);
            for(j=i+1;j<Ngr;++j)
                searcher.search(excluded_groups[i], excluded_groups[j], nlist);
        }
        searcher.search(excluded_groups[Ngr-1], not_excluded, nlist);
    }
    // Compute neighbours inside not excluded group
    searcher.search(not_excluded, nlist);
}
*/

///////////////////////////////////
inline float LJ_en_kernel(float sig, float eps, float r){
    float tmp = sig/r;
    tmp = tmp*tmp; // This gets (s/r)^2
    tmp = tmp*tmp*tmp; // This gets (s/r)^6
    return 4.0*eps*(tmp*tmp-tmp);
}

inline float Coulomb_en_kernel(float q1, float q2, float r){
    return ONE_4PI_EPS0*q1*q2/r;
}

// Apply combination rule
void Simulation::apply_comb_rule_GROMOS(int p1, int p2, float& sig, float& eps){
    float c6,c12;
    //GROMACS rule 1. C6-C12 provided, not sigma-epsilon!!!
    c6  = sqrt(sigma(p1)*sigma(p2));
    c12 = sqrt(epsilon(p1)*epsilon(p2));
    // Convert to sigma and epsilon
    c6 !=0 ? sig = pow(c12/c6,1.0/6.0) : sig = 0.0;
    c12!=0 ? eps = 0.25*c6*c6/c12 : eps = 0.0;
}

void Simulation::apply_comb_rule_AMBER(int p1, int p2, float& sig, float& eps){
    //GROMACS rule 2
    sig = 0.5*(sigma(p1)+sigma(p2));
    eps = sqrt(epsilon(p1)*epsilon(p2));
}

void Simulation::apply_comb_rule_OPLS(int p1, int p2, float& sig, float& eps){
    //GROMACS rule 3
    sig = sqrt(sigma(p1)*sigma(p2));
    eps = sqrt(epsilon(p1)*epsilon(p2));
}

// nb energy of the pair of atoms
Energy_components Simulation::non_bond_energy(int a1, int a2, Frame& target){
    std::vector<Eigen::Vector2i> nl(1);
    nl[0](0) = a1;
    nl[0](1) = a2;
    return non_bond_energy(nl,target);
}

// Interaction energy of two selections
Energy_components Simulation::non_bond_energy(Selection& sel1, Selection& sel2, float cutoff, bool is_periodic){
    std::vector<Eigen::Vector2i> nl;
    Grid_searcher(cutoff, sel1, sel2, nl, true, is_periodic);
    return non_bond_energy(nl,sel1.get_system()->Frame_data(sel1.get_frame()));
}

/*
// Total nb energy for global nlist
Energy_components Simulation::non_bond_energy(){
    return non_bond_energy(nlist);
}
*/

// Non-bond energy for given nlist
Energy_components Simulation::non_bond_energy(const std::vector<Eigen::Vector2i>& nlist,
                                              Frame& target_frame,
                                              std::vector<float>* dist_vec){
    int i,j,n,p1,p2;
    float scale_q, scale_lj;
    float r;
    float sig,eps; //actual LJ params after applying comb. rule
    float e1, e2;
    bool is_14;

    // Separate energies to report
    Energy_components en;

    n = nlist.size();
    for(i=0;i<n;++i){
        p1 = nlist[i](0); //We have absolute indexes in nlist!
        p2 = nlist[i](1);
        // If this is bonded pair skip
        if(bond_dist(p1,p2)<3) continue;
        // Check for 1-4 scaling
        if(bond_dist(p1,p2)==3){
            // 1-4 pair
            is_14 = true;
            scale_q = scale14_coulomb;
            scale_lj = scale14_lj;
        } else {
            // normal pair
            is_14 = false;
            scale_q = 1.0;
            scale_lj = 1.0;
        }

        // Get distance
        if(!dist_vec){
            r = (target_frame.coord[p2]-target_frame.coord[p1]).norm();
        } else {
            r = (*dist_vec)[i];
        }

        //cout << "+++ " << p1 << " " << p2 << endl;
        //cout << "+++ " << all.XYZ(p1).transpose() << " " << all.XYZ(p2).transpose() << endl;
        //cout << "+++ " << r << endl;
        //cout << "+++ " << all.XYZ(p2) << " " << all.XYZ(p1) << endl;

        // Compute Coulomb interaction
        e1 = scale_q*Coulomb_en_kernel(charge(p1), charge(p2), r);
        if(is_14)
            en.q_14 += e1;
        else
            en.q_sr += e1;

        // Call needed combination rule
        comb_rule_function(p1,p2,sig,eps);

        // Compute LJ interaction
        e2 = scale_lj*LJ_en_kernel(sig,eps,r);
        if(is_14)
            en.lj_14 += e2;
        else
            en.lj_sr += e2;

        en.total += e1+e2; // Add terms to total energy
    }

    return en;
}

/*
void Simulation::non_bond_force(){
    int i,j,n,p1,p2;
    float scale_q, scale_lj;
    float en = 0.0;
    float inv_r;
    float sig,eps;
    float e1, e2;
    Vector3f dist_vec;
    bool is_14;

    n = nlist.size();
    for(i=0;i<n;++i){
        p1 = nlist[i](0); //We have absolute indexes in nlist!
        p2 = nlist[i](1);
        // If this is bond pair skip
        if(bond_dist(p1,p2)<3) continue;
        // Check for 1-4 scaling
        if(bond_dist(p1,p2)==3){
            is_14 = true;
            scale_q = scale14_coulomb;
            scale_lj = scale14_lj;
        } else {
            is_14 = false;
            scale_q = 1.0;
            scale_lj = 1.0;
        }

        // Get distance vector
        dist_vec = all.XYZ(p2)-all.XYZ(p1);
        // Get inverse distance
        inv_r = 1.0/dist_vec.norm();

        // Compute Coulomb interaction
        e1 = scale_q*ONE_4PI_EPS0*charge(p1)*charge(p2) * inv_r*inv_r*inv_r;
        force.col(p1) -= dist_vec*e1;
        force.col(p2) += dist_vec*e1;

        // Compute LJ interaction
        apply_comb_rule(p1,p2,sig,eps);

        sig = sig * inv_r; //sigma/r
        sig = sig*sig; //(sigma/r)^2
        sig = sig*sig*sig; // (sigma/r)^6

        e2 = scale_lj*24.0*eps*inv_r*inv_r*(sig-2.0*sig*sig);
        force.col(p1) += dist_vec*e2;
        force.col(p2) -= dist_vec*e2;
    }
}

void Simulation::compute_force(){
    force.fill(0.0);
    non_bond_force();
}

void Simulation::verlet_step(){
    int i;
    float t = params.time_step;

    for(i=0;i<Natoms;++i){
        all.XYZ(i,0) += velocity.col(i)*t
                        + 0.5*t*t* force.col(i)*inv_mass(i);
        velocity.col(i) += 0.5*t* force.col(i)*inv_mass(i); //v_half
    }

    compute_force(); // Compute with new coordinates

    for(i=0;i<Natoms;++i)
        velocity.col(i) += 0.5*t* force.col(i)*inv_mass(i); // Second step
}

*/
