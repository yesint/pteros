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
#include <unordered_set>
#include <unordered_map>
#include <Eigen/Core>
#include <boost/function.hpp>



namespace pteros {

/**
  */

struct Force_field {
    std::vector<Eigen::Vector2i> charge_groups;
    std::vector<std::unordered_set<int>> exclusions;
    Eigen::MatrixXf LJ_C6, LJ_C12; // Matrices for normal (not excluded, not 1-4) LJ interactions
    std::vector<Eigen::Vector2f> LJ14_interactions; // [C6,C12]
    std::unordered_map<int,int> LJ14_pairs; // (a*LJ14_interactions.size()+b) --> LJ14_interaction
    float fudgeQQ; // 1-3 charge scaling factor
    bool ready;

    Force_field(){
        ready = false;
    }

    void clear(){
        charge_groups.clear();
        exclusions.clear();
        LJ_C6.fill(0.0);
        LJ_C12.fill(0.0);
        LJ14_interactions.clear();
        LJ14_pairs.clear();
        fudgeQQ = 0.0;
        ready = false;
    }
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
