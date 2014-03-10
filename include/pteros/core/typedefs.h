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

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <Eigen/Dense>

namespace pteros {

typedef Eigen::Ref<Eigen::Vector3f> Vector3f_ref;
typedef Eigen::Ref<Eigen::Matrix3f> Matrix3f_ref;
typedef Eigen::Ref<Eigen::VectorXf> VectorXf_ref;
typedef Eigen::Ref<Eigen::MatrixXf> MatrixXf_ref;
typedef Eigen::Ref<Eigen::Vector3i> Vector3i_ref;

typedef Eigen::Ref<const Eigen::Vector3f> Vector3f_cref;
typedef Eigen::Ref<const Eigen::Matrix3f> Matrix3f_cref;
typedef Eigen::Ref<const Eigen::VectorXf> VectorXf_cref;
typedef Eigen::Ref<const Eigen::MatrixXf> MatrixXf_cref;
typedef Eigen::Ref<const Eigen::Vector3i> Vector3i_cref;

}

#endif
