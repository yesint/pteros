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

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <Eigen/Dense>

namespace pteros {

using Vector3f_ref = Eigen::Ref<Eigen::Vector3f> ;
using Matrix3f_ref = Eigen::Ref<Eigen::Matrix3f> ;
using VectorXf_ref = Eigen::Ref<Eigen::VectorXf> ;
using MatrixXf_ref = Eigen::Ref<Eigen::MatrixXf> ;
using Vector3i_ref = Eigen::Ref<Eigen::Vector3i> ;

using Vector3f_const_ref = const Eigen::Ref<const Eigen::Vector3f>& ;
using Matrix3f_const_ref = const Eigen::Ref<const Eigen::Matrix3f>& ;
using VectorXf_const_ref = const Eigen::Ref<const Eigen::VectorXf>& ;
using MatrixXf_const_ref = const Eigen::Ref<const Eigen::MatrixXf>& ;
using Vector3i_const_ref = const Eigen::Ref<const Eigen::Vector3i>& ;

}

#endif
