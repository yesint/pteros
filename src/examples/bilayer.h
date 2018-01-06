/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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



#ifndef BILAYER_H
#define BILAYER_H

#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/distance_search.h"

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
    std::shared_ptr<Selection> spot1_ptr, spot2_ptr;
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

/*
float point_in_membrane(Eigen::Vector3f& point, Selection& head_markers, float d,
                        const Eigen::Vector3i& pbc_dims = Eigen::Vector3i::Ones());
*/
}

#endif

