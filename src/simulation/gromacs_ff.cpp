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
#include "pteros/core/grid_search.h"
#include "pteros/core/pteros_error.h"
#include <Eigen/Core>
#include <fstream>
#include <sstream>

using namespace std;
using namespace Eigen;

namespace pteros {
// Gromacs moleculetype
struct Moleculetype {
    vector<Atom> atoms;
    vector<Vector2i> bonds;
    vector<Vector3i> angles;
    vector<Vector4i,Eigen::aligned_allocator<Eigen::Vector4i> > dihedrals;
    vector<Vector2i> pairs;
};

// Gromacs atom type
struct Atom_type {
    std::string name;
    std::string type;
    float sigma, epsilon;
};

void Simulation::read_gromacs_top(std::string fname){
    ifstream f(fname.c_str());
    if(!f) throw Pteros_error("Failed to read top file "+fname);

    cout << "Reading processed GROMACS topology file '" << fname << "'..." << endl;

    Selection all(*system,"all");
    cout << all.size() << endl;

    stringstream ss;
    string line, dum, mode, mtype;
    float v;
    int comb;

    // Initialize bond distance matrix
    bond_dist.resize(Natoms,Natoms);
    bond_dist.fill(4); // Separated by more then 3 bonds by default

    // Stuctures for assigning atom types
    std::vector<Atom_type> atom_types;
    std::map<std::string,int> type_map;

    // Moleculetypes
    map<string,Moleculetype> moleculetypes;
    // Molecules
    map<string,int> molecules;

    // Prepare bonded arrays
    constr.clear();
    angles.clear();
    dihedrals.clear();

    while(getline(f,line)){
        if(line.size()<=1) continue; //Skip empty lines
        ss.clear(); ss.str(line);
        ss >> dum;
        if(dum[0]==';') continue; //Skip comments
        // If we are here, then this is an informative line
        if(dum[0]=='['){
            // This is a directive, parse it
            ss >> mode;
            cout << "Reading " << mode << "... " << endl;
            continue;
        }

        // Depending on the mode parse line
        ss.clear(); ss.str(line);
        if(mode=="defaults"){
            // Reade fudge factors and combination rule
            ss >> dum >> comb >> dum >> scale14_lj >> scale14_coulomb;
            switch(comb){
                case 1:                    
                    comb_rule_function = boost::bind(&Simulation::apply_comb_rule_GROMOS,this,_1,_2,_3,_4);
                    break;
                case 2:                    
                    comb_rule_function = boost::bind(&Simulation::apply_comb_rule_AMBER,this,_1,_2,_3,_4);
                    break;
                case 3:                    
                    comb_rule_function = boost::bind(&Simulation::apply_comb_rule_OPLS,this,_1,_2,_3,_4);
                    break;
                default:
                    throw Pteros_error("Unsupported combination rule!");
            }
        } else if(mode=="atomtypes"){
            Atom_type t;
            ss >> t.name >> t.type >> v >> v >> dum >> t.sigma >> t.epsilon;
            // Depending of the comb rule we read either c6 and c12 or sigma and epsilon!
            // For comb_rule==1 C6-C12 is read, beware!
            atom_types.push_back(t);
            type_map[t.name] = atom_types.size()-1;
        } else if(mode=="bondtypes"){
        } else if(mode=="angletypes"){
        } else if(mode=="dihedraltypes"){
        } else if(mode=="moleculetype"){
            ss >> mtype;
            moleculetypes[mtype] = Moleculetype();
            cout << "\t" << mtype << endl;
        } else if(mode=="atoms"){
            Atom at; // Temporary atom
            string type_str;
            float q,m;
            int i;
            ss >> i >> type_str >> dum >> dum >> dum >> dum >> q >> m;
            at.charge = q;
            at.mass = m;
            at.type = type_map[type_str];
            moleculetypes[mtype].atoms.push_back(at);
        } else if(mode=="bonds"){
            int i,j;
            ss >> i >> j;
            // Add to bond list
            moleculetypes[mtype].bonds.push_back(Vector2i(i-1,j-1));
        } else if(mode=="constraints"){
            int i,j;
            ss >> i >> j;
            // Add to bond list
            moleculetypes[mtype].bonds.push_back(Vector2i(i-1,j-1));
        } else if(mode=="pairs"){
            // 1-4 pairs for which modified non-bond forces should be applied
            int i,j;
            ss >> i >> j;
            moleculetypes[mtype].pairs.push_back(Vector2i(i-1,j-1));
        } else if(mode=="angles"){
            int i,j,k;
            ss >> i >> j >> k;
            moleculetypes[mtype].angles.push_back(Vector3i(i-1,j-1,k-1));
        } else if(mode=="atomtypes"){
        } else if(mode=="dihedrals"){
            int i,j,k,l;
            ss >> i >> j >> k >> l;
            moleculetypes[mtype].dihedrals.push_back(Vector4i(i,j,k,l));
        } else if(mode=="system"){
        } else if(mode=="molecules"){
            // Reading molecules
            string mol_name;
            int n;
            ss >> mol_name >> n;
            molecules[mol_name] = n;
            cout << "\t" << mol_name << " " << n << endl;
        }
    }

    f.close();

    cout << "Assigning moleculetypes to the system..." << endl;
    int last_atom = 0;
    int ind, at1, at2, at3, at4;
    // Assign moleculetypes to molecules and map everything to the pteros::system
    for(map<string,int>::iterator it = molecules.begin();it!=molecules.end();it++){
        int n = it->second;
        mtype = it->first;
        cout << "\tMolecule: " << mtype << " Num: " << n << endl;
        // For each molecule of this type do
        for(int m=0;m<n;++m){
            // Assign atoms
            for(int i=0; i<moleculetypes[mtype].atoms.size(); ++i){                
                ind = i + last_atom + m*moleculetypes[mtype].atoms.size();                
                all.Mass(ind) = moleculetypes[mtype].atoms[i].mass;
                all.Charge(ind) = moleculetypes[mtype].atoms[i].charge;
                all.Type(ind) = moleculetypes[mtype].atoms[i].type;
            }
            // Assign 1-4 pairs
            for(int i=0; i<moleculetypes[mtype].pairs.size(); ++i){
                at1 = moleculetypes[mtype].pairs[i](0) + last_atom + m*moleculetypes[mtype].atoms.size();
                at2 = moleculetypes[mtype].pairs[i](1) + last_atom + m*moleculetypes[mtype].atoms.size();
                // Add to bond_dist matrix
                bond_dist(at1,at2) = 3;
                bond_dist(at2,at1) = 3;
            }
            // Assign bonds
            for(int i=0; i<moleculetypes[mtype].bonds.size(); ++i){
                at1 = moleculetypes[mtype].bonds[i](0) + last_atom + m*moleculetypes[mtype].atoms.size();
                at2 = moleculetypes[mtype].bonds[i](1) + last_atom + m*moleculetypes[mtype].atoms.size();
                // Add to bond_dist matrix
                bond_dist(at1,at2) = 1;
                bond_dist(at2,at1) = 1;
                // Add to the list of constraints
                constr.push_back(Vector2i(at1,at2));
            }
            // Assign angles
            for(int i=0; i<moleculetypes[mtype].angles.size(); ++i){
                at1 = moleculetypes[mtype].angles[i](0) + last_atom + m*moleculetypes[mtype].atoms.size();
                at2 = moleculetypes[mtype].angles[i](1) + last_atom + m*moleculetypes[mtype].atoms.size();
                at3 = moleculetypes[mtype].angles[i](2) + last_atom + m*moleculetypes[mtype].atoms.size();
                // Add to bond_dist matrix
                bond_dist(at1,at3) = 2;
                bond_dist(at3,at1) = 2;
                // Add to the list of angles
                angles.push_back(Vector3i(at1,at2,at3));
            }
            // Assign dihedrals
            for(int i=0; i<moleculetypes[mtype].dihedrals.size(); ++i){
                at1= moleculetypes[mtype].dihedrals[i](0)+last_atom + m*moleculetypes[mtype].atoms.size();
                at2= moleculetypes[mtype].dihedrals[i](1)+last_atom + m*moleculetypes[mtype].atoms.size();
                at3= moleculetypes[mtype].dihedrals[i](2)+last_atom + m*moleculetypes[mtype].atoms.size();
                at4= moleculetypes[mtype].dihedrals[i](3)+last_atom + m*moleculetypes[mtype].atoms.size();
                // Add to the list of dihedrals
                dihedrals.push_back(Vector4i(at1,at2,at3,at4));
            }
        }
        // End of molecules of this type
        last_atom += n*moleculetypes[mtype].atoms.size();
    }

    cout << "Assigning non-bond parameters by atomtypes..." << endl;
    // Fill fast internal arrays
    for(int i=0;i<Natoms;++i){
        sigma(i) = atom_types[all.Type(i)].sigma;
        epsilon(i) = atom_types[all.Type(i)].epsilon;
        inv_mass(i) = 1.0/(all.Mass(i));
        charge(i) = all.Charge(i);
    }

    //cout << sigma.transpose() << endl;

    cout << "Topology is ready!" << endl << endl;
}

}
