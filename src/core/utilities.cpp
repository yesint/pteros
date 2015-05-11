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

#include "pteros/core/utilities.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

float pteros::angle_between_vectors(Vector3f_const_ref vec1, Vector3f_const_ref vec2)
{
    float ang = vec1.dot(vec2)/vec1.norm()/vec2.norm();
    if(ang>1.0) ang = 1.0;
    if(ang<-1.0) ang = -1.0;
    return acos(ang);
}


Vector3f project_vector(Vector3f_const_ref vec1, Vector3f_const_ref vec2)
{
    return (vec1.dot(vec2)/vec2.dot(vec2))*vec2;
}


float rad_to_deg(float rad)
{
    return rad * 57.295779513082320876798154814105;
}
