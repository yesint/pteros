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

#ifndef ATOM_H
#define ATOM_H

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
    /// Atom name (CA, O, N, etc.)
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
    /// @}

    /// @name Force-field related fields
    /// @{
    /// Atom mass
    float mass;
    /// Atom charge
    float charge;
    /// Atom type code. -1 means unknown.
    int type;
    /// Atom type name - textual representation of the atom type
    std::string type_name;
    /// @}

    Atom(){
        resid = -1;
        name = "";
        chain = ' ';
        resname = "";
        tag = "";
        occupancy = 0;
        beta = 0;
        resindex = -1;
        mass = 0;
        charge = 0;
        type = -1;
        type_name = "";        
    }
};

}

#endif /* ATOM_H */
