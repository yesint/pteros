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


struct Species {
    Species(): total_area(0.0), total_volume(0.0) {}

    Selection sel;
    int num_residues;
    double total_area;
    double total_volume;
    // Mapping from pid to local indexes
    map<int,int> pid_to_ind;
    // Mapping from local indexes to pid
    map<int,int> ind_to_pid;
};


struct Face {
    Face(int at1, int at2, double a): pids(at1,at2),area(a) {}
    Vector2i pids;
    double area;
};


void compute_voronoi_3d(const vector<Selection>& species_sel){
    using namespace voro;

    // Make selections for species
    vector<Species> species(species_sel.size());
    for(int i=0;i<species_sel.size();++i){
        species[i].sel = species_sel[i];
        species[i].num_residues = species[i].sel.get_resindex(true).size();
    }

    // Initialize stats
    vector<Face> area_elements;

    // Initialize container
    auto m = species[0].sel.box().get_matrix();
    container_periodic con(m.col(0)(0),
                           m.col(1)(0),m.col(1)(1),
                           m.col(2)(0),m.col(2)(1),m.col(2)(2),
                           10,10,10,8);

    // Add particles from all species to container and set mapping
    map<int,int> pid_to_species;
    int pid = 0;
    for(int j=0;j<species.size();++j){
        auto& sp = species[j];
        for(int i=0;i<sp.sel.size();++i){
            con.put(pid, sp.sel.x(i), sp.sel.y(i), sp.sel.z(i));
            sp.pid_to_ind[pid]=i;
            sp.ind_to_pid[i]=pid;
            pid_to_species[pid] = j;
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
        int cur_sp_ind = pid_to_species[cur_pid];

        // Areas
        vector<double> face_areas;
        c.face_areas(face_areas);
        vector<int> neib;
        c.neighbors(neib); // Get neigbour pids

        // Neighbors and faces arrays are consistent to each other
        for(int i=0;i<neib.size();++i){
            // Get resindex of neighbour
            int neib_pid = neib[i];
            auto& neib_sp_ind = pid_to_species[neib_pid];

            if(cur_sp_ind != neib_sp_ind){
                // This is a face between different species
                area_elements.emplace_back(cur_pid,neib_pid,face_areas[i]);
                species[cur_sp_ind].total_area += face_areas[i];
            }
        }
        // Add to volume of current species
        species[cur_sp_ind].total_volume += c.volume();

    } while(clo.inc());

    fmt::print("There are {} inter-species faces\n",area_elements.size());

    // Compute areas of all species-species interfaces
    MatrixXd interf(species.size(),species.size());
    MatrixXd interf_per_res(species.size(),species.size());
    interf.fill(0.0);
    interf_per_res.fill(0.0);
    for(const auto& el: area_elements){
        int sp1 = pid_to_species[el.pids[0]];
        int sp2 = pid_to_species[el.pids[1]];
        interf(sp1,sp2) += el.area;
        interf(sp2,sp1) += el.area;
    }
    // Diagonal element are species areas
    for(int i=0;i<species.size();++i) interf(i,i)=species[i].total_area;

    for(int i=0;i<species.size();++i){
        interf_per_res.col(i) = interf.col(i)/float(species[i].num_residues);
    }

    // Print it
    fmt::print("Species interfaces areas:\n{}\n",interf);
    fmt::print("Species interfaces areas per residue:\n{}\n",interf_per_res);

}

} // namespace pteros
