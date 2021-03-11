/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2021, Semen Yesylevskyy
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


#include "gromacs_utils.h"

using namespace pteros;
using namespace std;

namespace pteros {

void init_gmx_box(matrix box)
{
    box[0][0] = 0.0; box[0][1] = 0.0; box[0][2] = 0.0;
    box[1][0] = 0.0; box[1][1] = 0.0; box[1][2] = 0.0;
    box[2][0] = 0.0; box[2][1] = 0.0; box[2][2] = 0.0;
}


void gmx_box_to_pteros(const matrix m, PeriodicBox &b)
{
    Eigen::Matrix3f matr;
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            matr(i,j) = m[j][i];
    b.set_matrix(matr);
}

void pteros_box_to_gmx(const PeriodicBox &box, matrix m)
{
    Eigen::Matrix3f matr = box.get_matrix();
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            m[i][j] = matr(j,i);
}

}
