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


#ifndef BILAYER_H
#define BILAYER_H

#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/grid_search.h"

namespace pteros {

struct Bilayer_point_info {
    /// To which monolayer it belongs
    int monolayer;
    /// Distances from marker surfaces
    float surf_dist1, surf_dist2;
    /// Distance from central plane
    float center_dist;
    /// Vector of bilayer normal
    Eigen::Vector3f normal;
    /// Position of bilayer center
    Eigen::Vector3f center;
    /// Bilayer thickness
    float thickness;
    /// Projections of the point to both surfaces
    Eigen::Vector3f proj1, proj2;
    /// Pointers to projection spots (marker atoms in both monolayers used to compute normal)
    boost::shared_ptr<Selection> spot1_ptr, spot2_ptr;
    void print();
};

class Bilayer {
public:
    Bilayer(){};
    Bilayer(Selection& sel, std::string head_marker_atom, float d = 2.0);
    void create(Selection& sel, std::string head_marker_atom, float d = 2.0);
    /// Given a point in space computes various properties of bilayer in the vicinity of this point
    Bilayer_point_info point_info(Eigen::Vector3f& point);
protected:
    Selection* bilayer_ptr;
    Selection mono1,mono2; // Monolayers
    std::vector<Selection> surf; // Marker surfaces
    int spot_size;
};

}

#endif
