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

struct Face {
    Face(int at1, int at2, double a): atoms(at1,at2),area(a) {}
    Vector2i atoms;
    double area;
};


void compute_voronoi_3d(const Selection& sel, const string& solvent_sel){
    using namespace voro;


    // Initialize stats
    vector<Face> area_elements;

    LOG()->info("aaaa");

    auto m = sel.box().get_matrix();

    container_periodic con(m.col(0)(0),
                           m.col(1)(0),m.col(1)(1),
                           m.col(2)(0),m.col(2)(1),m.col(2)(2),
                           10,10,10,8);

    LOG()->info("bbb");

    // Add particles to container
    for(int i=0;i<sel.size();++i){
        con.put(i, sel.x(i), sel.y(i), sel.z(i));
    }

    // Compute cells
    c_loop_all_periodic clo(con);
    voronoicell_neighbor c;
    if(! clo.start()) throw PterosError("No particles!");
    do {
        con.compute_cell(c,clo);

        int cur_at = clo.pid(); // Current atom
        int cur_resind = sel.resindex(cur_at); // Current resindex

        // Areas
        vector<double> face_areas;
        c.face_areas(face_areas);
        vector<int> neib;
        c.neighbors(neib); // Get neigbour atoms local indexes

        // Neighbors and faces arrays are consistent to each other
        for(int i=0;i<neib.size();++i){
            int r = sel.resindex(neib[i]);
            if(r!=cur_resind && !(sel.resname(cur_at)=="SOL" && sel.resname(neib[i])=="SOL")){
                // This is a face between the residues
                area_elements.emplace_back(cur_at,neib[i],face_areas[i]);
            }
        }
    } while(clo.inc());

    fmt::print("There are {} inter-residue faces\n",area_elements.size());

    //con.draw_cells_gnuplot("cells.gnu");
    //con.draw_particles("p.dat");
    for(auto& el: area_elements){
        fmt::print("{} {} {}\n",sel.resname(el.atoms[0]),sel.resname(el.atoms[1]),el.area);
    }
}

}
