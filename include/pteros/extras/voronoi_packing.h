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
#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>

namespace pteros {

/*
class VoronoiPacking {
public:
    VoronoiPacking();
    // Adds new group and returns its id
    int add_group(const std::string& name, const Selection& sel);
    // All interface areas. Column i gives areas of group i with all groups
    Eigen::MatrixXd all_interface_areas();
    // Return a vector of interface areas between this group and other groups
    Eigen::VectorXd group_interface_areas(int gr);
    Eigen::VectorXd group_interface_areas(const std::string& gr);
    // Interface area between given groups
    double pair_interface_area(const std::string& gr1,const std::string& gr2);
    double pair_interface_area(int gr1,int gr2);
    // Group volume
    double group_volume(int gr);
    double group_volume(const std::string& gr);
    // Group area
    double group_area(int gr);
    double group_area(const std::string& gr);
    // Average inside the group by residues
    // Return a vector of average per residue type interface areas between this and other groups
    Eigen::VectorXd per_resname_interface_areas(int gr, const std::string& resname);
    // Per-atom stats for given resname
    Eigen::VectorXd per_resname_interface_areas(int gr, const std::string& resname);
private:
};
*/

void compute_voronoi_3d(const std::vector<Selection>& groups_sel);

} // namespace pteros




