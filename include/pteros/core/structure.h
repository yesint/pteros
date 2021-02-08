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
#include <Eigen/Core>
#include "pteros/core/atom.h"
#include "pteros/core/state.h"

namespace pteros {

class Structure {
public:
    Structure();

    size_t size() const {return _atoms.size();}
    size_t num_atoms() const {return _atoms.size();}
    size_t num_bonds() const {return _bonds.size();}
    size_t num_molecules() const {return _molecules.size();}
    bool has_bonds() const {return !_bonds.empty();}
    bool has_molecules() const {return !_molecules.empty();}

    void compute_bonds(const State& state, float cutoff=0.15);
    void compute_molecules();

    inline Atom& atom(size_t ind) { return _atoms[ind]; }

    std::vector<Atom>& atoms() {return _atoms;}
    std::vector<Eigen::Vector2i>& bonds() {return _bonds;}
    std::vector<Eigen::Vector2i>& molecules() {return _molecules;}

    std::vector<int> bonded_atoms(int at) const;
    std::vector<int> find_bonds(int at) const;
    int find_molecule(int at) const;

    void clear();

private:
    std::vector<Atom> _atoms;
    std::vector<Eigen::Vector2i> _bonds;
    std::vector<Eigen::Vector2i> _molecules;
}; // Structure


} // pteros
