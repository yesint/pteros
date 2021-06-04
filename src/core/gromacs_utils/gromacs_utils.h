/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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


#pragma once

#include "pteros/core/typedefs.h"
#include "pteros/core/periodic_box.h"
#include <Eigen/Core>

#ifdef USE_GROMACS
#include "gromacs/math/vectypes.h"
#else
#include "xdrfile.h"
#endif

namespace pteros {

void init_gmx_box(matrix box);
void gmx_box_to_pteros(const matrix m, PeriodicBox& b);
void pteros_box_to_gmx(const PeriodicBox& b, matrix m);

}
