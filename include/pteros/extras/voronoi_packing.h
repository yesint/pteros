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
#include <Eigen/Core>
//#include "voro++.hh"

namespace pteros {

void compute_voronoi_3d(const std::vector<Selection>& groups_sel);

struct PackingGroup {
    PackingGroup(): total_area(0.0), total_volume(0.0) {} 

    Selection sel;
    int num_residues;
    double total_area;
    double total_volume;
    // Mapping from pid to local indexes
    std::map<int,int> pid_to_ind;
    // Mapping from local indexes to pid
    std::map<int,int> ind_to_pid;

    void compute_averages(int num_frames);
};


struct InterGroupFace {
    InterGroupFace(int at1, int at2, double a): pids(at1,at2),area(a) {}
    Eigen::Vector2i pids;
    double area;
};


class Voronoi3D {
public:
    Voronoi3D(){}
    Voronoi3D(const std::vector<Selection>& groups_sel);
    void create(const std::vector<Selection>& groups_sel);
    void compute();
    void compute_averages();
    void write_stats(const std::string& fname) const;
    PackingGroup& get_group(int i) { return groups[i]; }
    int num_groups() const { return groups.size(); }
    Eigen::MatrixXd& get_interface_areas() { return interface_areas; }
    int get_num_frames() const { return num_frames; }
private:
    int num_frames;
    std::map<int, int> pid_to_groups;
    std::vector<PackingGroup> groups;
    // Output matrix of areas between groups
    Eigen::MatrixXd interface_areas;
};

} // namespace pteros




