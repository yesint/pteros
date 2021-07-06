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


#include "pteros/core/utilities.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include <Eigen/Core>
#include "voro++.hh"
#include <string>

using namespace std;
using namespace pteros;
using namespace Eigen;

namespace pteros {


struct PackingGroup {
    PackingGroup(): total_area(0.0), total_volume(0.0) {}

    Selection sel;
    int num_residues;
    double total_area;
    double total_volume;
    // Mapping from pid to local indexes
    map<int,int> pid_to_ind;
    // Mapping from local indexes to pid
    map<int,int> ind_to_pid;
};


struct InterGroupFace {
    InterGroupFace(int at1, int at2, double a): pids(at1,at2),area(a) {}
    Vector2i pids;
    double area;
};


void compute_voronoi_3d(const vector<Selection>& groups_sel){
    using namespace voro;

    // Make selections for species
    vector<PackingGroup> groups(groups_sel.size());
    for(int i=0;i<groups_sel.size();++i){
        groups[i].sel = groups_sel[i];
        groups[i].num_residues = groups[i].sel.get_resindex(true).size();
    }

    // Initialize stats
    vector<InterGroupFace> area_elements;

    // Initialize container
    auto m = groups[0].sel.box().get_matrix();
    container_periodic con(m.col(0)(0),
                           m.col(1)(0),m.col(1)(1),
                           m.col(2)(0),m.col(2)(1),m.col(2)(2),
                           10,10,10,8);

    // Add particles from all species to container and set mapping
    map<int,int> pid_to_groups;
    int pid = 0;
    for(int j=0;j<groups.size();++j){
        auto& sp = groups[j];
        for(int i=0;i<sp.sel.size();++i){
            con.put(pid, sp.sel.x(i), sp.sel.y(i), sp.sel.z(i));
            sp.pid_to_ind[pid]=i;
            sp.ind_to_pid[i]=pid;
            pid_to_groups[pid] = j;
            ++pid;
        }
    }

    fmt::print("Computing cells...");

    // Compute cells
    c_loop_all_periodic clo(con);
    voronoicell_neighbor c;
    if(! clo.start()) throw PterosError("No particles!");
    do {
        con.compute_cell(c,clo);

        int cur_pid = clo.pid(); // Current pid
        int cur_gr_ind = pid_to_groups[cur_pid];

        // Areas
        vector<double> face_areas;
        c.face_areas(face_areas);
        vector<int> neib;
        c.neighbors(neib); // Get neigbour pids

        // Neighbors and faces arrays are consistent to each other
        for(int i=0;i<neib.size();++i){
            // Get resindex of neighbour
            int neib_pid = neib[i];
            auto& neib_sp_ind = pid_to_groups[neib_pid];

            if(cur_gr_ind != neib_sp_ind){
                // This is a face between different species
                area_elements.emplace_back(cur_pid,neib_pid,face_areas[i]);
                groups[cur_gr_ind].total_area += face_areas[i];
            }
        }
        // Add to volume of current species
        groups[cur_gr_ind].total_volume += c.volume();

    } while(clo.inc());

    fmt::print("There are {} inter-group faces\n",area_elements.size());

    // Compute areas of all species-species interfaces
    MatrixXd interf(groups.size(),groups.size());
    MatrixXd interf_per_res(groups.size(),groups.size());
    interf.fill(0.0);
    interf_per_res.fill(0.0);
    for(const auto& el: area_elements){
        int sp1 = pid_to_groups[el.pids[0]];
        int sp2 = pid_to_groups[el.pids[1]];
        interf(sp1,sp2) += el.area;
        interf(sp2,sp1) += el.area;
    }
    // Diagonal element are species areas
    for(int i=0;i<groups.size();++i) interf(i,i)=groups[i].total_area;

    for(int i=0;i<groups.size();++i){
        interf_per_res.col(i) = interf.col(i)/float(groups[i].num_residues);
    }

    // Print it
    fmt::print("Group interfaces areas:\n{}\n",interf);
    fmt::print("Group interfaces areas per residue:\n{}\n",interf_per_res);

}

} // namespace pteros
