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


#include "pteros/core/system.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "pteros/core/force_field.h"
#include "pteros/core/utilities.h"
#include <cmath>
#include <functional>
#include "pteros/core/logging.h"


using namespace std;
using namespace pteros;
using namespace Eigen;

Vector3f get_shift_coefs(int alpha, float r1, float rc){
    Vector3f res;
    res(0) = -(( (alpha+4)*rc - (alpha+1)*r1 )/( pow(rc,alpha+2)*pow(rc-r1,2) ));
    res(1) = ( (alpha+3)*rc - (alpha+1)*r1 )/( pow(rc,alpha+2)*pow(rc-r1,3) );
    res(2) = 1.0/pow(rc,alpha) - (res(0)/3.0)*pow(rc-r1,3) - (res(1)/4.0)*pow(rc-r1,4);
    return res;
}


// Plain LJ kernel
float LJ_en_kernel(float C6, float C12, float r, const ForceField& ff){
    float r_inv = 1.0/r;
    float tmp = r_inv*r_inv; // (1/r)^2
    tmp = tmp*tmp*tmp; // (1/r)^6
    return C12*tmp*tmp-C6*tmp;
}

// Cutoff LJ kernel
float LJ_en_kernel_cutoff(float C6, float C12, float r, const ForceField& ff){
    if(r>ff.rvdw) return 0.0;
    float r_inv = 1.0/r;
    float tmp = r_inv*r_inv; // (1/r)^2
    tmp = tmp*tmp*tmp; // (1/r)^6
    return C12*tmp*tmp-C6*tmp;
}

// Shifted LJ kernel
float LJ_en_kernel_shifted(float C6, float C12, float r, const ForceField& ff){
    if(r>ff.rvdw) return 0.0;

    float val12 = pow(r,-12)
            -(ff.shift_12(0)/3.0)*pow(r-ff.rvdw_switch,3)
            -(ff.shift_12(1)/4.0)*pow(r-ff.rvdw_switch,4)
            -ff.shift_12(2);
    float val6 = pow(r,-6)
            -(ff.shift_6(0)/3.0)*pow(r-ff.rvdw_switch,3)
            -(ff.shift_6(1)/4.0)*pow(r-ff.rvdw_switch,4)
            -ff.shift_6(2);

    return C12*val12 - C6*val6;
}


#define ONE_4PI_EPS0      138.935456

// Plane Coulomb kernel
float Coulomb_en_kernel(float q1, float q2, float r, const ForceField& ff){
    return ff.coulomb_prefactor*q1*q2/r;
}

// Cutoff Coulomb kernel
float Coulomb_en_kernel_cutoff(float q1, float q2, float r, const ForceField& ff){
    if(r>ff.rcoulomb) return 0.0;
    return ff.coulomb_prefactor*q1*q2/r;
}


// Reaction field Coulomb kernel
float Coulomb_en_kernel_rf(float q1, float q2, float r, const ForceField& ff){
    return ff.coulomb_prefactor*q1*q2*(1.0/r + ff.k_rf*r*r - ff.c_rf);
}

// Shifted Coulomb kernel
float Coulomb_en_kernel_shifted(float q1, float q2, float r, const ForceField& ff){
    return ff.coulomb_prefactor*q1*q2*( 1.0/r
                                     -(ff.shift_1(0)/3.0)*pow(r-ff.rcoulomb_switch,3)
                                     -(ff.shift_1(1)/4.0)*pow(r-ff.rcoulomb_switch,4)
                                     -ff.shift_1(2)
                                     );
}


#define LOWER(s) str_to_lower_copy(s)

void ForceField::setup_kernels(){    
    LOG()->debug("Coulomb type: {}",coulomb_type);
    LOG()->debug("Coulomb modifier: {}",coulomb_modifier);
    LOG()->debug("VdW type: {}",vdw_type);
    LOG()->debug("VdW modifier: {}",vdw_modifier);

    // Set Coulomb prefactor
    coulomb_prefactor = ONE_4PI_EPS0 / epsilon_r;

    // Set Coulomb kernel
    if(LOWER(coulomb_type)=="reaction-field"){
        // In case of reaction field precompute constanst
        if(epsilon_rf){
            k_rf = (1.0/(rcoulomb*rcoulomb*rcoulomb))
                    * (epsilon_rf-epsilon_r) / (2.0*epsilon_rf+epsilon_r);
        } else {
            // for epsilon_rf = 0 (which means inf)
            k_rf = 0.5/(rcoulomb*rcoulomb*rcoulomb);
        }
        c_rf = (1.0/rcoulomb) + k_rf*rcoulomb*rcoulomb;

        // Set coulomb kernel pointer
        coulomb_kernel_ptr = &Coulomb_en_kernel_rf;
        LOG()->debug("\tCoulomb kernel: reaction_field");

    } else if( ( LOWER(coulomb_type)=="cut-off"
                 && LOWER(coulomb_modifier)== "potential-shift"
               )
              || LOWER(coulomb_type)=="shift"
              || LOWER(coulomb_type)=="pme") {
        // Compute shift constants for power 1
        shift_1 = get_shift_coefs(1,rcoulomb_switch,rcoulomb);

        coulomb_kernel_ptr = &Coulomb_en_kernel_shifted;
        LOG()->debug("\tCoulomb kernel: shifted");

    } else if(LOWER(coulomb_type)==LOWER("cut-off")) {
        // In other cases set plain Coulomb interaction
        coulomb_kernel_ptr = &Coulomb_en_kernel_cutoff;
        LOG()->debug("\tCoulomb kernel: cutoff");
    } else {
        coulomb_kernel_ptr = &Coulomb_en_kernel;
        LOG()->debug("\tCoulomb kernel: plain");
    }

    // Set LJ kernel
    if(LOWER(vdw_type)== "shift"){
        // Compute shift constants for powers 6 and 12
        shift_6 = get_shift_coefs(6,rvdw_switch,rvdw);
        shift_12 = get_shift_coefs(12,rvdw_switch,rvdw);

        LJ_kernel_ptr = &LJ_en_kernel_shifted;
        LOG()->debug("\tLJ kernel: shifted");

    } else if(LOWER(vdw_type)== "cut-off") {
        LJ_kernel_ptr = &LJ_en_kernel_cutoff;
        LOG()->debug("\tLJ kernel: cutoff");

    } else {
        LJ_kernel_ptr = &LJ_en_kernel;
        LOG()->debug("\tLJ kernel: plain");
    }
}


Vector2f ForceField::pair_energy(int at1, int at2, float r, float q1, float q2, int type1, int type2)
{
    float c6,c12;
    // indexes have to be in increasing order
    if(at1>at2){
        std::swap(at1,at2);
        std::swap(q1,q2);
        std::swap(type1,type2);
    }
    // Check if the pair is excluded
    if(exclusions[at1].count(at2)) return {0,0};
    // Check if this pair is 1-4 pair
    auto it = LJ14_pairs.find(at1*natoms+at2);
    if(it==std::end(LJ14_pairs)){
        // normal pair
        c6 = LJ_C6(type1,type2);
        c12 = LJ_C12(type1,type2);
        return {coulomb_kernel_ptr(q1,q2,r,*this), LJ_kernel_ptr(c6,c12,r,*this)};
    } else {
        // 1-4 pair
        c6 = LJ14_interactions[it->second](0);
        c12 = LJ14_interactions[it->second](1);
        return {coulomb_kernel_ptr(q1,q2,r,*this)*fudgeQQ, LJ_kernel_ptr(c6,c12,r,*this)};
    }
}

float ForceField::get_cutoff(){
    return std::min(rcoulomb,rvdw);
}

ForceField::ForceField():  ready(false) {}

ForceField::ForceField(const ForceField &other){
    exclusions = other.exclusions;
    molecules = other.molecules;
    bonds = other.bonds;

    LJ_C6 = other.LJ_C6;
    LJ_C12 = other.LJ_C12;
    LJ14_interactions = other.LJ14_interactions;
    LJ14_pairs = other.LJ14_pairs;
    fudgeQQ = other.fudgeQQ;
    rcoulomb = other.rcoulomb;
    epsilon_r = other.epsilon_r;
    epsilon_rf = other.epsilon_rf;
    rcoulomb_switch = other.rcoulomb_switch;
    rvdw_switch = other.rvdw_switch;
    rvdw = other.rvdw;
    coulomb_type = other.coulomb_type;
    coulomb_modifier = other.coulomb_modifier;
    vdw_type = other.vdw_type;
    vdw_modifier = other.vdw_modifier;

    ready = other.ready;

    if(ready) setup_kernels();
}

ForceField &ForceField::operator=(ForceField other){    
    exclusions = other.exclusions;
    molecules = other.molecules;
    bonds = other.bonds;

    LJ_C6 = other.LJ_C6;
    LJ_C12 = other.LJ_C12;
    LJ14_interactions = other.LJ14_interactions;
    LJ14_pairs = other.LJ14_pairs;
    fudgeQQ = other.fudgeQQ;
    rcoulomb = other.rcoulomb;
    epsilon_r = other.epsilon_r;
    epsilon_rf = other.epsilon_rf;
    rcoulomb_switch = other.rcoulomb_switch;
    rvdw_switch = other.rvdw_switch;
    rvdw = other.rvdw;
    coulomb_type = other.coulomb_type;
    coulomb_modifier = other.coulomb_modifier;
    vdw_type = other.vdw_type;
    vdw_modifier = other.vdw_modifier;

    ready = other.ready;

    if(ready) setup_kernels();

    return *this;
}

void ForceField::clear(){    
    exclusions.clear();
    LJ_C6.fill(0.0);
    LJ_C12.fill(0.0);
    LJ14_interactions.clear();
    LJ14_pairs.clear();
    fudgeQQ = 0.0;
    molecules.clear();
    bonds.clear();

    ready = false;
}
