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
#include "pteros/core/typedefs.h"

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
 Matrix3f b = box.get_matrix();
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
    Periodic_box(Matrix3f_const_ref box);

    /// Constructor from vectors and angles
    Periodic_box(Vector3f_const_ref vectors, Vector3f_const_ref angles);

    /// Copy constructor
    Periodic_box& operator=(Periodic_box other){
        modify(other._box);
        return *this;
    }

    /// Modify the box
    void modify(Matrix3f_const_ref box);       

    /// Get i-th box vector
    Eigen::Vector3f get_vector(int i){ return _box.col(i); }

    /// Get stored matrix of box vectors
    Eigen::Matrix3f get_matrix() const {return _box;}

    /// Get stored inverted matrix of box vectors
    Eigen::Matrix3f get_inv_matrix() const {return _box_inv;}

    /// Convert point from lab coordinates to box coordinates
    Eigen::Vector3f lab_to_box(Vector3f_const_ref point) const
    { return _box_inv.colwise().normalized()*point; }

    /// Return the transformation from lab coordinates to box coordinates
    Eigen::Matrix3f lab_to_box_transform() const {return _box_inv.colwise().normalized();}

    /// Convert point from box coordinates to lab coordinates
    Eigen::Vector3f box_to_lab(Vector3f_const_ref point) const
    { return _box.colwise().normalized()*point; }

    /// Return the transformation from box coordinates to lab coordinates
    Eigen::Matrix3f box_to_lab_transform() const {return _box.colwise().normalized();}

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
    float distance(Vector3f_const_ref point1,
                   Vector3f_const_ref point2,         
                   Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// The same as distance but returns squared distance
    float distance_squared(Vector3f_const_ref point1,
                   Vector3f_const_ref point2,                   
                   Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// Wrap point to the box for given set of dimensions
    /// Origin of the box coordinates is assumed to be {0,0,0}.
    void wrap_point(Vector3f_ref point,
                    Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// Determine if the point is inside the box
    /// Origin of the box coordinates is assumed to be {0,0,0}.
    bool in_box(Vector3f_const_ref point);

    /// Finds a periodic image of point, which is closest in space to target and returns it    
    Eigen::Vector3f get_closest_image(Vector3f_const_ref point,
                                      Vector3f_const_ref target,
                                      Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// Computes shortest vector from point1 to point2 between their closest images
    Eigen::Vector3f shortest_vector(Vector3f_const_ref point1,
                                    Vector3f_const_ref point2,
                                    Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// Returns box volume
    float volume();

    /// Read box from CRYST string in PDB format. Overwrites current box!
    void read_pdb_box(const char *line);

    /// Write box as CRYST string in PDB format
    std::string write_pdb_box() const;

    /// Returns representation of the box as direction vectors and angles
    void to_vectors_angles(Vector3f_ref vectors, Vector3f_ref angles) const;

    /// Creates box from vectors and angles. Overwrites current box!
    void from_vectors_angles(Vector3f_const_ref vectors, Vector3f_const_ref angles);


private:
    Eigen::Matrix3f _box;
    Eigen::Matrix3f _box_inv;
    Eigen::Vector3f _extents;
    bool _is_triclinic;
    bool _is_periodic;
};

}

#endif /* ATOM_H */
