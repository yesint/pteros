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

#include "pteros/core/selection.h"
#include "pteros/analysis/frame_info.h"

namespace pteros {

class JumpRemover {
public:    
    JumpRemover();
    void add_atoms(const Selection &sel);
    void set_pbc(Array3i_const_ref pbc);
    /// -1 means do not unwrap, 0 means auto find distance
    void set_unwrap_dist(float d);
    void set_pbc_atom(int ind);

    // Remove jumps
    void remove_jumps(System& system);

private:    
    // Indexes for removing jumps
    std::vector<int> no_jump_ind;
    // Running reference coordinates for removing jumps
    Eigen::MatrixXf no_jump_ref;
    // Dimensions to consider
    Eigen::Array3i dims;
    // Starting distance for unwrapping. -1 means no unwrapping (default)
    float unwrap_d;
    // Leading index for unwrapping
    int pbc_atom;

    bool initialized;
};

} // namespace
