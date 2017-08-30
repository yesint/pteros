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

#ifndef UTILITIES_H
#define UTILITIES_H

#include "pteros/core/typedefs.h"

namespace pteros {

    float angle_between_vectors(Vector3f_const_ref vec1, Vector3f_const_ref vec2);

    Eigen::Vector3f project_vector(Vector3f_const_ref vec1, Vector3f_const_ref vec2);

    float rad_to_deg(float rad);

    constexpr long double operator"" _deg ( long double deg )
    {
        return deg*3.141592/180;
    }

/*
    float distance_to_vector(Vector3f_const_ref point,
                             Vector3f_const_ref vec_origin,
                             Vector3f_const_ref vec_dir);

    void project_to_vector(Vector3f_const_ref point,
                             Vector3f_const_ref vec_origin,
                             Vector3f_const_ref vec_dir);

    float distance_to_plane(Vector3f_const_ref point,
                            Vector3f_const_ref plane_origin,
                            Vector3f_const_ref plane_normal);
*/

}

#endif
