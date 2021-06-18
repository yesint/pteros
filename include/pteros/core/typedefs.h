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

using Array3i_const_ref = const Eigen::Ref<const Eigen::Array3i>& ;
using ArrayXi_const_ref = const Eigen::Ref<const Eigen::ArrayXi>& ;

static const Eigen::Array3i fullPBC = Eigen::Array3i::Ones();
static const Eigen::Array3i noPBC   = Eigen::Array3i::Zero();

}




