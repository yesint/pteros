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

#ifndef FORCE_FIELD_H
#define FORCE_FIELD_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <Eigen/Core>

namespace pteros {

/**
  Force field parameters of the system.
  MD packages usually separate topology and force filed i.e.
  interaction parameters and chemical identity and connectivity (bonds, angles, atom types, etc).
  This doesn't make much sense since atom types are part of particular force field and at the
  same time part of chemical identity, etc.
  In Pteros the force field is all information, which is needed to compute interactomic
  energies and forces in addition to what is stored already in Atom.

  Atom already contains the following information: mass, charge, type and type_name.
  All the rest is stored in this class.

*/

class Force_field {
public:
    /// Charge groups. Currently not used and could be deleted later.
    std::vector<Eigen::Vector2i> charge_groups;
    /// Exclusions.
    /// The storage order is: (atom)-->[i1,i2,i3...in]
    /// which means that all interactions (atom:i1), (atom:i2) ... (atom:in) are excluded
    /// If atom has no exclusions the set is empty
    std::vector<std::unordered_set<int> > exclusions;
    /// Matrices of normal (not excluded, not 1-4) LJ interactions.
    /// The size of the matrix == the number of distinct LJ types.
    /// The types themselves are integers and are stored in Atom.type
    /// Matrices are symmetric.
    Eigen::MatrixXf LJ_C6, LJ_C12;
    /// The list of distinct types of LJ14 interactions in the format [C6,C12]
    std::vector<Eigen::Vector2f> LJ14_interactions;
    /// The list of LJ14 pairs. Each (a,b) pair is encoded as (a*LJ14_interactions.size()+b)
    /// There is a mapping:
    /// (a*LJ14_interactions.size()+b) --> LJ14_interaction
    std::unordered_map<int,int> LJ14_pairs;
    /// Scaling factor of 1-4 Coulomb interactions
    float fudgeQQ;

    float rcoulomb, epsilon_r, epsilon_rf, rcoulomb_switch, rvdw_switch, rvdw;
    std::string coulomb_type, coulomb_modifier, vdw_type, vdw_modifier;

    /// Topology-related infromation
    std::vector<Eigen::Vector2i> bonds;

    /// Is the force field properly set up?
    bool ready;

    /// Pointer to chosen coulomb kernel
    std::function<float(float,float,float)> coulomb_kernel_ptr;

    /// Pointer to chosen VDW kernel
    std::function<float(float,float,float)> LJ_kernel_ptr;

    // Aux constants to be precomputed by set_kernels()
    float coulomb_prefactor, k_rf, c_rf;
    // potential shift constants
    Eigen::Vector3f shift_1, shift_6, shift_12;

    /// Constructor
    Force_field();

    /// Copy constructor
    Force_field(const Force_field& other);

    /// Assignment operator
    Force_field& operator=(Force_field other);

    // Clear ff
    void clear();

    // Setup coulomb and VDW kernel pointers
    void setup_kernels();

    float LJ_en_kernel(float C6, float C12, float r);
    float LJ_en_kernel_cutoff(float C6, float C12, float r);
    float LJ_en_kernel_shifted(float C6, float C12, float r);
    float Coulomb_en_kernel(float q1, float q2, float r);
    float Coulomb_en_kernel_rf(float q1, float q2, float r);
    float Coulomb_en_kernel_shifted(float q1, float q2, float r);
    float Coulomb_en_kernel_cutoff(float q1, float q2, float r);
};


}

#endif /* FORCE_FIELD_H */
