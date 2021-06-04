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

#include <string>
#include <Eigen/Core>
#include "pteros/core/typedefs.h"

namespace pteros {

/** @brief Class encapsulating all operations with arbitrary triclinic periodic boxes
 This class stores the periodic box itself and also contains pre-computed matrices
 for converting laboratory coordinates to the coordinates in triclinic basis and vice versa.
 Extents of the periodic box are also precomputed and stored internally. 
 */
class PeriodicBox {
public:

    /// Default constructor
    PeriodicBox();

    /// Constructor from matrix
    PeriodicBox(Matrix3f_const_ref m);

    /// Constructor from vector lengths and angles
    PeriodicBox(Vector3f_const_ref vectors, Vector3f_const_ref angles);

    /// Copy constructor
    PeriodicBox(const PeriodicBox& other);

    /// Assignment operator
    PeriodicBox& operator=(const PeriodicBox& other);

    // Get box element
    float get_element(int i, int j) const;

    // Set box element
    void set_element(int i, int j, float val);

    /// Get i-th box vector
    Eigen::Vector3f get_vector(int i) const;

    /// Set i-th box vector
    void set_vector(Vector3f_const_ref vec, int i);

    /// Get stored matrix of box vectors
    Eigen::Matrix3f get_matrix() const;

    /// Modify the box from 3x3 matrix
    void set_matrix(Matrix3f_const_ref box);

    /// Scale box vectors by specified factors.
    /// Causes recomputing internal data.
    void scale_vectors(Vector3f_const_ref scale);

    /// Get stored inverted matrix of box vectors
    Eigen::Matrix3f get_inv_matrix() const;

    /// Convert point from lab coordinates to box coordinates
    Eigen::Vector3f lab_to_box(Vector3f_const_ref point) const;

    /// Return the transformation from lab coordinates to box coordinates
    Eigen::Matrix3f lab_to_box_transform() const;

    /// Convert point from box coordinates to lab coordinates
    Eigen::Vector3f box_to_lab(Vector3f_const_ref point) const;

    /// Return the transformation from box coordinates to lab coordinates
    Eigen::Matrix3f box_to_lab_transform() const;

    /// Return i-th extent of the box
    float extent(int i) const;

    /// Return the vector of box extents
    Eigen::Vector3f extents() const;

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
                   Array3i_const_ref pbc = fullPBC) const;

    /// The same as distance but returns squared distance
    float distance_squared(Vector3f_const_ref point1,
                   Vector3f_const_ref point2,                   
                   Array3i_const_ref pbc = fullPBC) const;

    /// Wrap point to the box for given set of dimensions
    /// Origin of the box coordinates defaults to {0,0,0}.
    void wrap_point(Vector3f_ref point,
                    Array3i_const_ref pbc = fullPBC,
                    Vector3f_const_ref origin = Eigen::Vector3f::Zero()) const;

    /// Determine if the point is inside the box
    /// Origin of the box coordinates defaults to {0,0,0}.
    bool in_box(Vector3f_const_ref point, Vector3f_const_ref origin = Eigen::Vector3f::Zero()) const;

    /// Finds a periodic image of point, which is closest in space to target and returns it    
    Eigen::Vector3f closest_image(Vector3f_const_ref point,
                                      Vector3f_const_ref target,
                                      Array3i_const_ref pbc = fullPBC) const;

    /// Computes shortest vector from point1 to point2 between their closest images
    Eigen::Vector3f shortest_vector(Vector3f_const_ref point1,
                                    Vector3f_const_ref point2,
                                    Array3i_const_ref pbc = fullPBC) const;

    /// Returns box volume
    float volume();

    /// Read box from CRYST string in PDB format. Overwrites current box!
    void from_pdb_box(const char *line);

    /// Write box as CRYST string in PDB format
    std::string to_pdb_box() const;

    /// Returns representation of the box as vector lengths and angles
    void to_vectors_angles(Vector3f_ref vectors, Vector3f_ref angles) const;

    /// Creates box from vector length and angles. Overwrites current box!
    /// vectors = {a,b,c}
    /// angles = {a^c, b^c, a^b}
    void from_vectors_angles(Vector3f_const_ref vectors, Vector3f_const_ref angles);


private:
    Eigen::Matrix3f _box;
    Eigen::Matrix3f _box_inv;    
    bool _is_triclinic;
    bool _is_periodic;

    void recompute_internals();
};

}




