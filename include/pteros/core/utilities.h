/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2020, Semen Yesylevskyy
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

#include "pteros/core/typedefs.h"
#include "pteros/core/selection.h"
#include <vector>

namespace pteros {

    float angle_between_vectors(Vector3f_const_ref vec1, Vector3f_const_ref vec2);

    Eigen::Vector3f project_vector(Vector3f_const_ref vec1, Vector3f_const_ref vec2);

    float rad_to_deg(float ang);
    float deg_to_rad(float ang);

    constexpr long double operator"" _deg ( long double ang ) {
        return ang*3.141592/180.0;
    }

    constexpr long double operator"" _rad ( long double ang ) {
        return ang*180.0/3.141592;
    }


    std::string get_element_name(int elnum);

    int get_element_number(const std::string& name);

    float get_vdw_radius(int elnum, const std::string& name);

    void guess_element(const std::string& name, int& anum, float& mass);


    /// Returns rotation matrix given pivot, axis and angle in radians
    Eigen::Affine3f rotation_transform(Vector3f_const_ref pivot, Vector3f_const_ref axis, float angle);

    /// Get 1-letter protein code from 3-letter
    char resname_1char(const std::string& code);
    /// Get 3-letter protein code from 1-letter
    std::string resname_3char(char code);

    /// Simple histogram class
    class Histogram {
    public:
        Histogram(){}
        Histogram(float minval, float maxval, int n);
        void create(float minval, float maxval, int n);
        int get_bin(float v);
        void add(float v, float weight=1.0);
        void add(const std::vector<float> &v);
        void add_cylindrical(float r, float w, float sector, float cyl_h);
        void normalize(float norm=0);
        float value(int i) const;
        float position(int i) const;
        float delta() const {return d;}
        Eigen::VectorXd &values();
        Eigen::VectorXd& positions();
        int num_bins() const;
        void save_to_file(const std::string& fname, float x_shift=0);
    private:
        int nbins;
        float minv,maxv,d;
        Eigen::VectorXd val;
        Eigen::VectorXd pos;
        bool normalized;
    };

    /// Simple 2D histogram class
    class Histogram2D {
    public:
        Histogram2D(){}
        Histogram2D(float minval1, float maxval1, int n1, float minval2, float maxval2, int n2);
        void create(float minval1, float maxval1, int n1, float minval2, float maxval2, int n2);
        void add(float v1, float v2, float weight=1.0);
        void normalize(float norm=0);
        Eigen::Vector2f delta() const {return d;}

        void save_to_file(const std::string& fname);
    private:
        Eigen::Vector2i nbins;
        Eigen::Vector2f minv,maxv,d;
        Eigen::MatrixXd val;

        bool normalized;
    };

    void greeting(std::string tool_name="");

} // namespace



