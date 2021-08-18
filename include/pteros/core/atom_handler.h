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

#include "pteros/core/atom.h"
#include <Eigen/Core>

namespace pteros {


using AtomStateIterator = std::vector<Eigen::Vector3f>::iterator;
using AtomStateConstIterator = std::vector<Eigen::Vector3f>::const_iterator;
using AtomIterator = std::vector<Atom>::iterator;
using AtomConstIterator = std::vector<Atom>::const_iterator;

class Selection;
class System;

/// Auxilary type used to incapsulate the atom and its current coordinates
/// Used internally in Selection::operator[] and in iterator access to Selection.
/// Objects of this class should not be created by the user in normal situation.
class AtomHandler {
    friend class Selection;
    friend class System;
public:
    AtomHandler(){}
    AtomHandler(const Selection& sel, int i);
    AtomHandler(const System& sys, int i, int fr);
    void set(const Selection& sel, int i);
    void set(const System& sys, int i, int fr);    

    /// @name Inline accessors. Const and non-const versions.
    /// @{
    inline const int& index() const { return ind; }

    inline int& resid(){ return atom_ptr->resid; }
    inline const int& resid() const { return atom_ptr->resid; }

    inline std::string& name(){ return atom_ptr->name; }
    inline const std::string& name() const { return atom_ptr->name; }

    inline char& chain(){ return atom_ptr->chain; }
    inline const char& chain() const { return atom_ptr->chain; }

    inline std::string& resname(){ return atom_ptr->resname; }
    inline const std::string& resname() const { return atom_ptr->resname; }

    inline std::string& tag(){ return atom_ptr->tag; }
    inline const std::string& tag() const { return atom_ptr->tag; }

    inline float& occupancy(){ return atom_ptr->occupancy; }
    inline const float& occupancy() const { return atom_ptr->occupancy; }

    inline float& beta(){ return atom_ptr->beta; }
    inline const float& beta() const { return atom_ptr->beta; }

    inline int& resindex() { return atom_ptr->resindex; }
    inline const int& resindex() const { return atom_ptr->resindex; }

    inline float& mass(){ return atom_ptr->mass; }
    inline const float& mass() const { return atom_ptr->mass; }

    inline float& charge(){ return atom_ptr->charge; }
    inline const float& charge() const { return atom_ptr->charge; }

    inline int& type(){ return atom_ptr->type; }
    inline const int& type() const { return atom_ptr->type; }

    inline std::string& type_name(){ return atom_ptr->type_name; }
    inline const std::string& type_name() const { return atom_ptr->type_name; }

    inline float& x(){ return (*coord_ptr)(0); }
    inline const float& x() const { return (*coord_ptr)(0); }

    inline float& y(){ return (*coord_ptr)(1); }
    inline const float& y() const { return (*coord_ptr)(1); }

    inline float& z(){ return (*coord_ptr)(2); }
    inline const float& z() const { return (*coord_ptr)(2); }

    inline Eigen::Vector3f& xyz(){ return *coord_ptr; }
    inline const Eigen::Vector3f& xyz() const { return *coord_ptr; }

    /*
    inline Eigen::Vector3f& vel(){ return *v_ptr; }
    inline const Eigen::Vector3f& vel() const { return *v_ptr; }

    inline Eigen::Vector3f& force(){ return *f_ptr; }
    inline const Eigen::Vector3f& force() const { return *f_ptr; }
    */

    inline Atom& atom(){ return *atom_ptr; }
    inline const Atom& atom() const { return *atom_ptr; }

    inline int& atomic_number(){ return atom_ptr->atomic_number; }
    inline const int& atomic_number() const { return atom_ptr->atomic_number; }

    std::string element_name() const;

    float vdw() const;
    ///@}

private:
    AtomIterator atom_ptr;
    AtomStateIterator coord_ptr;
    // Atom don't store it's own index, so store it here
    int ind;

    // Move handler by n atoms. For internal usage.
    inline void advance(int n=1){
        coord_ptr+=n;
        atom_ptr+=n;
        ind+=n;
    }
};

} //namespace



