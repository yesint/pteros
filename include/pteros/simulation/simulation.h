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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "pteros/core/selection.h"
#include "pteros/core/grid_search.h"

// Class for simulation involving all-atom MD parameters
namespace pteros {

// User-supplied simulation parameters
class MD_params {
public:
    float nb_cutoff;
    bool is_periodic;
    float temperature;
    float time_step;

    MD_params(){
        nb_cutoff = 1.0;
        is_periodic = true;
        temperature = 300;
        time_step = 0.002;
    }
};

/// Components of the non-bond energy
struct Energy_components {
    /// Total energy
    float total;
    /// Lenard-Jones energy of 1-4 pairs
    float lj_14;
    /// Coloumb energy of 1-4 pairs
    float q_14;
    /// Short-range Lenard-Jones energy (within cut-off)
    float lj_sr;
    /// Short-range Coloumb energy (within cut-off)
    float q_sr;

    Energy_components(){
        total = 0.0;
        lj_14 = 0.0;
        q_14 = 0.0;
        lj_sr = 0.0;
        q_sr = 0.0;
    }
    /// Writes all energy components to string in the following order:
    /// total, lj_sr, lj_14, q_sr, q_14
    std::string to_str();
};

/// @brief Simulation class. This class is HIGHLY EXPERIMENTAL!
///
/// API is subject to change and the general design may become absolutely different.
/// Now only non_bond_energy() function is used by the stable part of Pteros.
class Simulation {
public:
    /// Default constructor
    Simulation();
    /// Constructs simulation from system and topology file
    Simulation(System& sys, std::string top_file);

    ~Simulation();

    /// Copy constructor
    Simulation(const Simulation& other);

    /// Assignment operator
    Simulation& operator=(Simulation other);

    /// Load structure and topology
    void create(System& sys, std::string top_file);

    /// Set working frame
    void set_frame(int fr){
        frame = fr;
    }


    /*
        /// Add excluded groups
        void add_excluded_groups(std::vector<std::string>& gr);
        /// Setup the system and check all parameters
        void setup();        
        /// Updates nlist
        void update_nlist();
        */
    /*
        /// Non-bond energy of the whole system
        Energy_components non_bond_energy();
*/
    /// Non-bond energy for externally given neighbour list
    /// If dist_vec is given then distances are not computed but taken from dist_vec
    Energy_components non_bond_energy(const std::vector<Eigen::Vector2i>& nlist,
                                      Frame& target_frame,
                                      std::vector<float>* dist_vec = NULL);
    /// Non-bond energy for externally given pair of atoms
    Energy_components non_bond_energy(int a1, int a2, Frame& target);
    /// Non-bond interaction energy of two selections
    Energy_components non_bond_energy(Selection& sel1, Selection& sel2,
                                      float cutoff = 1.0, bool is_periodic = false);
    /*
        void non_bond_force();
        /// Set working frame
        void set_frame(int fr){ all.set_frame(fr); }

        /// Computes force
        void compute_force();
        /// Do one time step of Verlet integration
        void verlet_step();

        /// Simulation parameters
        MD_params params;

        void get_force(Eigen::MatrixXf& f){
            f = force;
        }
*/
protected:
    /// Parent system
    System* system;
    /// Working frame
    int frame;

    /// Number of atoms in the system
    int Natoms;

    /// Velocities
    Eigen::MatrixXf velocity;
    /// Forces
    Eigen::MatrixXf force;
    /// Inverse masses (for speed)
    Eigen::VectorXf inv_mass;
    /// Lenard-Jones params for atoms
    Eigen::VectorXf sigma, epsilon;
    /// Charges
    Eigen::VectorXf charge;
    /*
        /// Neighbour list containing all non-bond interacting pairs in the system
        std::vector<Eigen::Vector2i> nlist;
*/
//    /// Grid searcher object
//    Grid_searcher searcher;

    /// Non-bond parameters
    /// Pointer to combination rule function
    boost::function<void(int,int,float&,float&)> comb_rule_function;

    /// 1-4 scaling factors
    float scale14_coulomb, scale14_lj;

    /// Bonded parameters
    /// The list of constraints (constrained bonds)
    std::vector<Eigen::Vector2i> constr;
    /// The list of angles
    std::vector<Eigen::Vector3i> angles;
    /// The list of digedrals
    // Use Eigen allocator to ensure correct alignment
    std::vector<Eigen::Vector4i,Eigen::aligned_allocator<Eigen::Vector4i> > dihedrals;

    /// Bond distance map.
    /// bond_dist[i,j] shows how many bonds separate atoms i and j.
    Eigen::MatrixXi bond_dist;

    /*
        /// Excluded groups
        std::vector<Selection> excluded_groups;
        /// Atoms outside of any excluded group
        Selection not_excluded;

        /// Setup should be verified before run
        bool setup_ok;
        */

    /// Topology reader
    void read_gromacs_top(std::string fname);

    /// Functions for different combination rules
    void apply_comb_rule_GROMOS(int p1, int p2, float& sig, float& eps);
    void apply_comb_rule_AMBER(int p1, int p2, float& sig, float& eps);
    void apply_comb_rule_OPLS(int p1, int p2, float& sig, float& eps);

};

} //namespace
#endif //SIMULATION_H
