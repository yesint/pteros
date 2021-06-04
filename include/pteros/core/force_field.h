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


#pragma once

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

class ForceField {
public:
    int natoms;
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

    /// The list of LJ14 pairs.
    /// The mapping is (a*Natoms+b)->index in LJ14_interactions
    /// a<=b! This is important
    std::unordered_map<int,int> LJ14_pairs;

    /// Scaling factor of 1-4 Coulomb interactions
    float fudgeQQ;

    float rcoulomb, epsilon_r, epsilon_rf, rcoulomb_switch, rvdw_switch, rvdw;
    std::string coulomb_type, coulomb_modifier, vdw_type, vdw_modifier;

    /// Bonds
    std::vector<Eigen::Vector2i> bonds;

    // Molecules
    std::vector<Eigen::Vector2i> molecules;

    /// Is the force field properly set up?
    bool ready;

    /// Pointer to chosen coulomb kernel
    float (*coulomb_kernel_ptr)(float,float,float,const ForceField&);

    /// Pointer to chosen VDW kernel
    float (*LJ_kernel_ptr)(float,float,float,const ForceField&);

    // Aux constants to be precomputed by set_kernels()
    float coulomb_prefactor, k_rf, c_rf;
    // potential shift constants
    Eigen::Vector3f shift_1, shift_6, shift_12;

    /// Constructor
    ForceField();

    /// Copy constructor
    ForceField(const ForceField& other);

    /// Assignment operator
    ForceField& operator=(ForceField other);

    // Clear ff
    void clear();

    // Setup coulomb and VDW kernel pointers
    void setup_kernels();

    // Computes energy of atom pair at given distance
    // Returns {Coulomb_en,LJ_en}
    Eigen::Vector2f pair_energy(int at1, int at2, float r, float q1, float q2, int type1, int type2);

    /// Returns "natural" cutoff (currenly min of rcoulomb and rvdw)
    float get_cutoff();
};

}




