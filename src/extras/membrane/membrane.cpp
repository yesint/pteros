/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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

#include "pteros/extras/membrane.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <ctime>

using namespace std;
using namespace pteros;
using namespace Eigen;



Membrane::Membrane(System *sys, const std::vector<Lipid_descr> &species): system(sys), lipid_species(species)
{
    log = std::make_shared<spdlog::logger>("membrane", Log::instance().console_sink);
    log->set_pattern(Log::instance().generic_pattern);

    Selection all_pos_sel(*system);

    // Creating selections
    for(auto& sp: lipid_species){
        vector<Selection> res;
        system->select(sp.whole).split_by_residue(res);
        for(auto& lip: res){
            auto mol = Lipid(lip,sp);
            all_pos_sel.append(mol.get_position_sel());
            lipids.push_back(mol);
            resindex_map[mol.get_position_sel().Resindex(0)] = lipids.size()-1;
        }
    }

    // Compute connectivity
    all_pos_sel.split_by_connectivity(2.0,leaflets,true);

    // Print statictics
    log->info("Number of lipids: {}",lipids.size());
    log->info("Number of leaflets: {}",leaflets.size());

}

Lipid::Lipid(const Selection &sel, const Lipid_descr &descr){
    name = descr.name;
    whole = sel;
    head_marker = whole(descr.head_marker);
    tail_marker = whole(descr.tail_marker);
    position_marker = whole(descr.position_marker);
}

