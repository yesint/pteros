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

#include "pteros/core/selection.h"
#include <Eigen/Core>

namespace pteros {

/** Implementation of the Gaussian Network Model (GNM) protein model.
  */
class GNM {
public:
    GNM(const Selection& sel, float cutoff=0.7);

    Eigen::VectorXf get_eigenvector(int i) const;
    Eigen::VectorXf get_B_factor() const;

    void write_eigenvectors(std::string fname, int v1, int v2);
    void compute_c_matrix();
    void compute_p_matrix();
    void write_c_matrix(std::string fname);
    void write_p_matrix(std::string fname);

    Eigen::MatrixXf get_subset_c_matrix(ArrayXi_const_ref subset) const;
    Eigen::MatrixXf get_c_matrix() const;
private:
    int N;
    Selection* sel;
    Eigen::MatrixXf eigenvectors;
    Eigen::VectorXf eigenvalues;
    Eigen::MatrixXf c,p;
    Eigen::VectorXf b; // B-factor
};

}




