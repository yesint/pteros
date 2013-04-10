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

#ifndef GNM_H_INCLUDED
#define GNM_H_INCLUDED

#include "pteros/core/selection.h"
#include <Eigen/Core>

namespace pteros {

/** Implementation of th Gaussian Network Model (GNM) protein model.
  */
class GNM {
    public:
        int N;
        Eigen::MatrixXf eigenvectors;
        Eigen::VectorXf eigenvalues;
        Eigen::MatrixXf c,p;

        GNM(Selection& sel, float cutoff);

        void compute(Selection& sel, float cutoff);
        void write_eigenvectors(std::string fname, int v1, int v2);
        void compute_c_matrix(bool normalize=false);
        void compute_p_matrix();
        void write_c_matrix(std::string fname);
        void write_p_matrix(std::string fname);
};

}
#endif // GNM_H_INCLUDED
