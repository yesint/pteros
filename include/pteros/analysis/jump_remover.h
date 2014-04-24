/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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

#ifndef JUMP_REMOVER_H
#define JUMP_REMOVER_H

#include "pteros/core/selection.h"
#include "pteros/analysis/frame_info.h"

#include <iostream>

namespace pteros {

class Jump_remover {
public:
    void add_no_jump_atoms(const Selection &sel);
    void remove_jumps(System& system, const Frame_info &info);
private:
    // Indexes for removing jumps
    std::vector<int> no_jump_ind;
    // Running reference coordinates for removing jumps
    Eigen::MatrixXf no_jump_ref;
};

} // namespace

#endif
