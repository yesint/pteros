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




#include "pteros/extras/solvate.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

namespace pteros {

std::map<std::string,int> solvate(System& solute, float d, std::string solvent_file, std::string custom_sel){
    System solvent;

    if(solvent_file==""){
        // Look for $GMXDATA environmental variable
        if (const char* env_gmx = std::getenv("GMXDATA")) {
            solvent_file = string(env_gmx)+"/top/spc216.gro";
        } else {
            throw PterosError("Can't find default gromacs solvent file spc216.gro!");
        }
    }

    LOG()->info("Loading solvent from '{}'...", solvent_file);
    solvent.load( solvent_file );

    if(solvent.box(0).is_triclinic())
        throw PterosError("Only rectangular solvent boxes are allowed!");

    // See how many solvent boxes should be used to cover solute

    // In arbitrary triclinic boxes find the maximal box coordinate
    Vector3f max_solute_coord = solute.box(0).box_to_lab( solute.box(0).extents() );
    Vector3f max_solvent_coord = solvent.box(0).extents();

    Vector3i nbox;
    for(int i=0; i<3; ++i) nbox(i) = int(ceil(max_solute_coord(i)/max_solvent_coord(i)));

    LOG()->info("Will use {} solvent boxes...", nbox.transpose());

    // Distribute solvent boxes
    {
        auto all = solvent.select_all();
        auto m = solvent.box(0).get_matrix();
        LOG()->info("Distributing solvent boxes...");
        solvent.distribute(all,nbox,m);
    }

    // Move min coords of solvent and solute to zero
    Vector3f solvent_min,solvent_max, solute_min, solute_max;
    auto solute_all = solute.select_all();
    auto solvent_all = solvent.select_all();

    solute_all.minmax(solute_min,solute_max);
    solvent_all.minmax(solvent_min,solvent_max);
    solvent_all.translate(-solvent_min);
    solute_all.translate(-solute_min);

    LOG()->info("Finding solvent atoms outside the solute box...");

    // Cut solvent atoms outside the solute box
    vector<int> bad;
    for(int i=0; i<solvent_all.size(); ++i){
        if( !solute.box(0).in_box(solvent_all.xyz(i)) ) bad.push_back(solvent_all.index(i));
    }

    // Select bad atoms
    Selection bad_sel(solvent,bad);
    // Select whole bad residues
    vector<Selection> bad_res;
    bad_sel.each_residue(bad_res);

    LOG()->info("Found {} solvent molecules outside the solute box...", bad_res.size());
    for(auto& sel: bad_res){
        sel.set_beta(-1000);
    }


    // Find last index of solute
    int last_solute_ind = solute.num_atoms()-1;

    // append good solvent to solute
    solute.append(solvent("beta > -1000"));

    // select overlapping water
    string s = fmt::format("by residue within {} pbc noself of index 0-{}", d, last_solute_ind);

    Selection sel(solute, s);

    LOG()->info("Found {} overlaping solvent atoms at cutoff={}", sel.size(),d);

    // Remove overlapping water
    sel.set_beta(-1000);

    // If we have custom selection use it
    if(custom_sel!=""){
        Selection sel(solute, custom_sel);
        LOG()->info("Removing atoms from custom selection '{}' ({} atoms)", custom_sel, sel.size());
        sel.set_beta(-1000);
    }

    // Translate back to initial box center

    solute.keep("beta > -1000");
    solute_all = solute();
    solute_all.translate(solute_min);

    // Report number of remaining solvent residues
    map<string,int> residues;
    int at=last_solute_ind+1;
    do {
        string resname = solute_all.resname(at);
        int resind = solute_all.resindex(at);

        // Find the end of this residue
        do {
            ++at;
        } while( at<solute_all.size() && solute_all.resindex(at) == resind);

        if(residues.count(resname)){
            // such resname is present
            ++residues[resname];
        } else {
            // new resname
            residues[resname] = 1;
        }

    } while(at<solute_all.size());

    LOG()->info("Number of solvent molecules added:");
    for(auto& it: residues){
        LOG()->info("\t{}: {}", it.first,it.second);
    }

    return residues;
}

} //namespace



