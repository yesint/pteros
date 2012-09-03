/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
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

/** Functions for working with periodic box entries in PDB files.
These functions are adapted from GROMACS 4.0 and modified heavily.
They are used to read and write arbitrary triclinic box from/to PDB fles.
*/


#ifndef PDB_BOX_H
#define PDB_BOX_H

#include <Eigen/Core>
#include <string>

namespace pteros {

/// Reads periodic box from string
void read_pdb_box(const char *line, Eigen::Matrix3f& box);
/// Writes periodic box to string
std::string write_pdb_box(const Eigen::Matrix3f& box);
/// Returns representation of the periodic box as vectors and angles
void box_to_vectors_angles(const Eigen::Matrix3f& box, Eigen::Vector3f& vectors, Eigen::Vector3f& angles);
}
#endif
