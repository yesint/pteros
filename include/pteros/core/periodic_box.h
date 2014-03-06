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

/** @brief Class encapsulating all operations with arbitrary triclinic periodic boxes
 This class stores the periodic box itself and also contains pre-computed matrices
 for converting laboratory coordinates to the coordinates in triclinic basis and vice versa.
 Extents of the periodic box are also precomputed and stored internally.
 All data in the class are read-only. The user can set the box by using the constructor
 or by calling modify(box), than all internal data would be precomputed.
 Individual components of the box can't be changed. The only way to change is to get the
 whole box, modify the component and set it back:
 \code
 Periodic_box box(some_data);
 Matrix3f b = box.get_box();
 b(1,2) *= 2.0; // Modify the box element
 box.modify(b); // Set new box. This will recompute internal matrices
 \endcode
 */
class Periodic_box {
public:

    /// Default constructor
    Periodic_box(){
        _is_periodic = false;
        _is_triclinic = false;
        _box.fill(0.0);
    }

    /// Constructor from other box
    Periodic_box(const Eigen::Matrix3f& box);

    /// Modify the box
    void modify(const Eigen::Matrix3f& box);

    /// Get stored periodic box
    Eigen::Matrix3f get_box() const {return _box;}

    /// Convert point from lab coordinates to box coordinates
    Eigen::Vector3f to_box(const Eigen::Vector3f& point) const { return _to_box*point; }

    /// Return the transformation matrix from lab coordinates to box coordinates
    const Eigen::Matrix3f& to_box_matrix() const {return _to_box;}

    /// Convert point from box coordinates to lab coordinates
    Eigen::Vector3f to_lab(const Eigen::Vector3f& point) const { return _to_lab*point; }

    /// Return the transformation matrix from box coordinates to lab coordinates
    const Eigen::Matrix3f& to_lab_matrix() const {return _to_lab;}

    /// Return i-th extent of the box
    float extent(int i) const {return _extents(i);}

    /// Return the vector of box extents
    const Eigen::Vector3f& extents() const {return _extents;}

    /// Is the box triclinic?
    bool is_triclinic() const {return _is_triclinic;}

    /// Is the box set? If false, the system is not periodic
    bool is_periodic() const {return _is_periodic;}

    /// Compute a periodic distance between two points in the box
    /// Periodicity is only accouted for given set of dimensions
    /// If you need a non-periodic distance over all dimensions it is more efficient
    /// to compute it directly as:
    ///\code
    /// float dist = (point2-point1).norm();
    ///\endcode
    float distance(const Eigen::Vector3f& point1,
                   const Eigen::Vector3f& point2,
                   bool do_wrapping = true,
                   const Eigen::Vector3i& periodic_dims = Eigen::Vector3i::Ones()) const;

    /// Wrap point to the box for given set of dimensions
    void wrap_point(Eigen::Vector3f& point,
                    const Eigen::Vector3i& dims_to_wrap = Eigen::Vector3i::Ones()) const;

    /// Finds a periodic image of point, which is closest in space to target and returns it
    /// This method wraps both point and targer to periodic box internally (this is usually what you want).
    /// If this is not needed set @param do_wrapping to false, but in this case make sure
    /// that wrapping is done manually before! Otherwise results would be incorrect.
    Eigen::Vector3f get_closest_image(const Eigen::Vector3f& point,
                                      const Eigen::Vector3f& target,
                                      bool do_wrapping = true,
                                      const Eigen::Vector3i& dims_to_wrap = Eigen::Vector3i::Ones()) const;

    /// Return box volume
    float volume();

    /// Read box from CRYST string in PDB format
    void read_pdb_box(const char *line);

    /// Write box as CRYST string in PDB format
    std::string write_pdb_box() const;

    /// Returns representation of the box as direction vectors and angles
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
