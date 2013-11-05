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

#ifndef PERIODIC_BOX_H
#define PERIODIC_BOX_H

#include <string>
#include <Eigen/Core>

namespace pteros {

class Periodic_box {
public:
    Periodic_box(const Eigen::Matrix3f& box);
    void modify(const Eigen::Matrix3f& box);

    Eigen::Matrix3f get_box() const {return _box;}

    Eigen::Vector3f to_box(const Eigen::Vector3f& point) const;
    const Eigen::Matrix3f& to_box_matrix() const {return _to_box;}

    Eigen::Vector3f to_lab(const Eigen::Vector3f& point) const;
    const Eigen::Matrix3f& to_lab_matrix() const {return _to_lab;}

    float extent(int i) const {return _extents(i);}
    const Eigen::Vector3f& extents() const {return _extents;}

    bool is_triclinic() const {return _is_triclinic;}
    bool is_periodic() const {return _is_periodic;}

    float distance(const Eigen::Vector3f& point1,
                   const Eigen::Vector3f& point2,
                   bool do_wrapping = true,
                   const Eigen::Vector3i& periodic_dims = Eigen::Vector3i::Ones()) const;

    void wrap_point(Eigen::Vector3f& point,
                    const Eigen::Vector3i& dims_to_wrap = Eigen::Vector3i::Ones()) const;

    Eigen::Vector3f get_closest_image(const Eigen::Vector3f& point,
                                      const Eigen::Vector3f& target,
                                      bool do_wrapping = true,
                                      const Eigen::Vector3i& dims_to_wrap = Eigen::Vector3i::Ones()) const;

    float volume();

    void read_pdb_box(const char *line);

    std::string write_pdb_box() const;

    void to_vectors_angles(Eigen::Vector3f& vectors, Eigen::Vector3f& angles) const;

private:
    Eigen::Matrix3f _box;
    Eigen::Matrix3f _to_box;
    Eigen::Matrix3f _to_lab;
    Eigen::Vector3f _extents;
    bool _is_triclinic;
    bool _is_periodic;
};

}

#endif /* ATOM_H */
