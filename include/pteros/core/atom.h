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

#include <string>

namespace pteros {
/**
* Class which represents a single atom.
* Coordinates are stored separately in the System.
*/
class Atom {
  public:
    /// @name General fields
    /// @{

    /// Residue ID (unique only inside given chain)
    int  resid;
    /// %Atom name (CA, O, N, etc.)
    std::string  name;
    /// Chain. Single letter (A,B,Z)
    char  chain;
    /// Residue name in 3-letters code (ALA, GLY, MET, etc.)
    std::string  resname;
    /// Arbitrary textual tag
    std::string  tag;
    /// Occupancy field
    float  occupancy;
    /// B-factor field
    float  beta;
    /// Unique residue ID.
    /// For multi-chain proteins this marks each residue in unique way
    /// regardless of the chain
    int  resindex;
    /// Element atomic number
    /// zero means particle, which is not real atom (for example dummy)
    int atomic_number;
    /// @}

    /// @name Force-field related fields
    /// @{

    /// %Atom mass
    float mass;
    /// %Atom charge
    float charge;
    /// %Atom type code. -1 means unknown.
    int type;
    /// %Atom type name - textual representation of the atom type
    std::string type_name;
    /// @}

    Atom():
        resid(-1),
        name(""),
        chain(' '),
        resname(""),
        tag(""),
        occupancy(0),
        beta(0),
        resindex(-1),
        mass(0),
        charge(0),
        type(-1),
        type_name(""),
        atomic_number(0)
    {}
};

}




