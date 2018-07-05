/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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

#include <Eigen/Core>
#include <vector>
#include <memory>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <pteros/core/selection.h>

void minmax_gpu(const std::vector<Eigen::Vector3f>& coords, const std::vector<int>& index);

class GPUSelection_impl;

class GPUSelection {
public:
    GPUSelection(const pteros::Selection& sel);
    Eigen::Vector3f center();
    virtual ~GPUSelection();
private:
    GPUSelection_impl* impl;
};

