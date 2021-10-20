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


#include "pteros/extras/voronoi_packing.h"
#include "pteros/core/utilities.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include <Eigen/Core>
#include "voro++.hh"
#include <string>
#include <fstream>
#include <thread>

using namespace std;
using namespace pteros;
using namespace Eigen;

namespace pteros {

Voronoi3D::Voronoi3D(const std::vector<Selection>& groups_sel)
{
    create(groups_sel);
}

void Voronoi3D::create(const std::vector<Selection> &groups_sel)
{
    num_frames = 0;
    // Make selections for groups and set mapping
    groups.resize(groups_sel.size());
    int pid = 0;
    for(int g=0;g<groups.size();++g){
        auto& gr = groups[g];
        gr.sel = groups_sel[g];
        gr.num_residues = gr.sel.get_resindex(true).size();
        for(int i=0;i<gr.sel.size();++i){
            gr.pid_to_ind[pid]=i;
            gr.ind_to_pid[i]=pid;
            pid_to_groups[pid] = g;
            ++pid;
        }
    }

    // Init interface areas matrix
    interface_areas.resize(groups.size(),groups.size());
    interface_areas.fill(0.0);
}

void Voronoi3D::compute(){    
    // Initialize container    
    auto m = groups[0].sel.box().get_matrix();
    voro::container_periodic con(m.col(0)(0),
                                 m.col(1)(0),m.col(1)(1),
                                 m.col(2)(0),m.col(2)(1),m.col(2)(2),
                                 10,10,10,8);

    // Update all selections for current frame and add particles to container    
    int pid = 0;
    for(auto& gr: groups){
        gr.sel.apply();
        for(int i=0;i<gr.sel.size();++i){
            con.put(pid, gr.sel.x(i), gr.sel.y(i), gr.sel.z(i));
            ++pid;
        }
    }
    
    // Compute cells
    voro::c_loop_all_periodic clo(con);
    voro::voronoicell_neighbor c;
    vector<InterGroupFace> area_elements;

    LOG()->info("Computing...");

    if(!clo.start()) throw PterosError("No particles!");
    do {
        con.compute_cell(c,clo.ijk,clo.q);

        int cur_pid = clo.pid(); // Current pid
        int cur_gr_ind = pid_to_groups[cur_pid];
        
        vector<double> face_areas;
        c.face_areas(face_areas); // Get areas
        vector<int> neib;
        c.neighbors(neib); // Get neigbour pids

        // Neighbors and faces arrays are consistent to each other
        for(int i=0;i<neib.size();++i){
            // Get resindex of neighbour
            int neib_pid = neib[i];
            auto& neib_gr_ind = pid_to_groups[neib_pid];

            if(cur_gr_ind != neib_gr_ind){
                // This is a face between different species
                area_elements.emplace_back(cur_pid,neib_pid,face_areas[i]);
                groups[cur_gr_ind].total_area += face_areas[i];
            }
        }
        // Add to volume of current group
        groups[cur_gr_ind].total_volume += c.volume();

    } while(clo.inc());    

    fmt::print("There are {} inter-group faces\n",area_elements.size());
    
    // Compute areas of all group-group interfaces
    for(const auto& el: area_elements){
        int gr1 = pid_to_groups[el.pids[0]];
        int gr2 = pid_to_groups[el.pids[1]];
        interface_areas(gr1,gr2) += el.area;
        interface_areas(gr2,gr1) += el.area;
    }
    // Diagonal element are total areas of groups for convenience
    for(int i=0;i<groups.size();++i) interface_areas(i,i) += groups[i].total_area;

    // Update number of frames consumed    
    ++num_frames;

    LOG()->info("Done");
}

void Voronoi3D::compute_averages()
{
    // Average all properties by number of consumed frames
    for(auto& gr: groups){
        gr.compute_averages(num_frames);
    }
    interface_areas /= float(num_frames);
}


void Voronoi3D::write_stats(const std::string& fname) const {
    // Output groups stats
    ofstream out(fname+"_groups.dat");
    fmt::print(out,"group\tarea\tvolume\ta_per_res\tv_per_res\n");
    for(int i=0; i<groups.size(); ++i){
        fmt::print(out,"{i}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n",i,
                   groups[i].total_area, groups[i].total_volume,
                   groups[i].total_area/float(groups[i].num_residues),
                   groups[i].total_volume/float(groups[i].num_residues)
                  );
    }
    out.close();

    // Output areas matrix
    out.open(fname+"_interfaces.dat");
    fmt::print(out,"{}",interface_areas);
    out.close();
}


void PackingGroup::compute_averages(int num_frames)
{
    total_area /= float(num_frames);
    total_volume /= float(num_frames);
}

} // namespace pteros
