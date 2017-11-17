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

    float rad_to_deg(float ang);
    float deg_to_rad(float ang);

    constexpr long double operator"" _deg ( long double ang )
    {
        return ang*3.141592/180.0;
    }

    constexpr long double operator"" _rad ( long double ang )
    {
        return ang*180.0/3.141592;
    }


    class Histogram {
    public:
        Histogram(float minval, float maxval, int n);
        void add(float v);
        void normalize();
        float value(int i) const;
        float position(int i) const;
        const Eigen::VectorXd& values() const;
        const Eigen::VectorXd& positions() const;
        int num_bins() const;
        void save_to_file(const std::string& fname);
    private:
        int nbins;
        float minv,maxv,d;
        Eigen::VectorXd val;
        Eigen::VectorXd pos;
        bool normalized;
    };





}

#endif
