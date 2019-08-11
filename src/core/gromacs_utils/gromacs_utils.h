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
void gmx_box_to_pteros(const matrix m, Periodic_box& b);
void pteros_box_to_gmx(const Periodic_box& b, matrix m);

}
