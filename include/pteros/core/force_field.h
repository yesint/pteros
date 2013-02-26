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

#ifndef FORCE_FIELD_H
#define FORCE_FIELD_H

#include <string>
#include <vector>
#include <Eigen/Core>
#include <boost/function.hpp>



namespace pteros {

// ==========================================================
//  UNIVERSAL things, which are the same for all force fields
// ==========================================================

/** Type of non-bond function
  Currently only NB_FUNC_LJ is supported
*/
enum NB_FUNC_TYPE {NB_FUNC_LJ, NB_FUNC_BUCKINGHAM};

/// Combination rules:
enum COMB_RULE_TYPE {
    /** rule type 1 in GROMACS.
        V == C6 = 4*epsilon*sigma^6
        W == C12 = 4*epsilon*sigma^12
        Vij = sqrt(Vi * Vj)
        Wij = sqrt(Wi * Wj)
    */
    COMB_RULE_GROMOS,
    /** rule type 2 in GROMACS
        V == sigma
        W == epsilon
        Vij = 0.5*(Vi + Vj)
        Wij = sqrt(Wi * Wj)
    */
    COMB_RULE_AMBER_CHARMM,
    /** rule type 3 in GROMACS
        V == sigma
        W == epsilon
        Vij = sqrt(Vi * Vj)
        Wij = sqrt(Wi * Wj)
    */
    COMB_RULE_OPLS
};

enum BOND_FUNC_TYPE {
    /** Normal harmonic bond (type 1)
        Vb(r) = 0.5*kb*(r-b)^2
        F(R) = kb*(r-b)*R/r
      */
    BOND_FUNC_BOND,
    /** 4-order bond used in GROMOS-96 ff (type 2)
        Vb(r) = (1/4)*kb*(r^2-b^2)^2
        F(R) = kb*(r^2-b^2)*R
      */
    BOND_FUNC_G96_BOND,
    /// Not currently supported
    BOND_FUNC_MORSE,
    BOND_FUNC_CUBIC,
    BOND_FUNC_CONNECTION,
    BOND_FUNC_HARMONIC_POT,
    BOND_FUNC_FENE,
    BOND_FUNC_TAB1,
    BOND_FUNC_TAB2,
    BOND_FUNC_RESTR
};


enum PAIRS_FUNC_TYPE {
    PAIRS_FUNC_1,
    PAIRS_FUNC_2
};


enum ANGLE_FUNC_TYPE {
    /** Harmonic angle, GROMACS type 1
      V_ijk(theta) = 0.5*k*(theta-theta0)^2
      */
    ANGLE_FUNC_ANGLE,
    /** Cosine-based angle, GROMACS type 2
      V_ijk(theta) = 0.5*k*(cos(theta)-cos(theta0))^2
      */
    ANGLE_FUNC_G96_ANGLE,
    ANGLE_FUNC_CROSS1,
    ANGLE_FUNC_CROSS2,
    /** Uray-Bradley potential
      V_ijk(theta) = 0.5*k*(theta-theta0)^2 + 0.5*kUB*(r_ik-r_ik0)^2
      */
    ANGLE_FUNC_UB,
    ANGLE_FUNC_QUARTIC,
    ANGLE_FUNC_TAB
};

enum DIHEDRAL_FUNC_TYPE {
    DIHEDRAL_FUNC_PROPER,
    DIHEDRAL_FUNC_IMPROPER,
    DIHEDRAL_FUNC_RB,
    DIHEDRAL_FUNC_PERIODIC_IMPROPER,
    DIHEDRAL_FUNC_FOURIER,
    DIHEDRAL_FUNC_TAB = 8,
    DIHEDRAL_FUNC_PROPER_MULTI = 9
};

enum CONSTRAINT_FUNC_TYPE {
    CONSTRAINT_FUNC_1,
    CONSTRAINT_FUNC_2
};


//------------------------------------------------

struct Interaction_func {
    std::string type;
    std::vector<float> params;
    boost::function<float(const std::vector<int>&, const std::vector<float>&)> energy;
    boost::function<std::vector<float>(const std::vector<int>&, const std::vector<float>&)> force;
};

struct Interaction {
    std::vector<int> ind;
    int func;
};

/**
  This is the force field of one particular system.
  It is read from the dump of GROMACS tpr file and contains
  preprocessed and prepared topology and interactions.
  Such approach is chosen because parsing GROMACS topologies is very difficult -
  one need to deal with all complex rules and subtleties, with parameter overrides
  in itp files, etc.

  Atom names, types, masses and charges are assigned directly to the System.atoms
  All other data are stored in this class.
  */

struct Force_field {
    std::vector<Eigen::Vector2i> charge_groups;
    std::vector<Interaction_func> interaction_funcs;
    std::vector<Interaction> interactions;
    /** Pair interactions matrix
        If element is 0 it is excluded
        If it is >0 corresponding interaction_func is used
      */
     Eigen::MatrixXi nb_interactions;
};


}

/*
//get type of atom i of molec type j
int A=mtop->moltype[j].atoms.atom[i].type;
//get type of atom l of molec type k
int B=mtop->moltype[k].atoms.atom[l].type;
t_functype interaction_type=mtop->ffparams.functype[B*atnr+A];
switch(interaction_type)
{
 case F_LJ:
  float c6=mtop->ffparams.iparams[B*atnr+A].lj.c6;
  float c12=mtop->ffparams.iparams[B*atnr+A].lj.c12;
  break;
 case F_MORSE:
  ....
}
 */

#endif /* FORCE_FIELD_H */
