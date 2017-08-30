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

#include "pteros/analysis/jump_remover.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;


Jump_remover::Jump_remover():
    dims(Eigen::Vector3i::Ones()),
    unwrap_d(-1.0),
    leading_index(0),
    initialized(false)
{ }

void Jump_remover::add_atoms(const Selection &sel)
{    
    int ind;
    int n = sel.get_system()->num_atoms();
    for(int i=0;i<sel.size();++i){
        ind = sel.Index(i);
        if(ind<0 || ind>=n) throw Pteros_error("Index {} for jump removal out of range (0:{})!",ind,n-1);
        no_jump_ind.push_back(ind);
    }

    // Only keep unique atoms!
    sort(no_jump_ind.begin(),no_jump_ind.end());
    auto it = unique(no_jump_ind.begin(), no_jump_ind.end());
    no_jump_ind.resize( it - no_jump_ind.begin() );
}

void Jump_remover::set_dimensions(Vector3i_const_ref dim)
{
    dims = dim;
    if(dims.sum()==0) LOG()->warn("No periodic dimensions, skipping jump removing.");
}

void Jump_remover::set_unwrap_dist(float d)
{
    unwrap_d = d;
}

void Jump_remover::set_leading_index(int ind)
{
    leading_index = ind;
}

void Jump_remover::remove_jumps(System& system){
    // Exit immediately if no atoms or no valid dimensions
    if(no_jump_ind.empty() || dims.sum()==0) return;
    // If not periodic also do nothing
    if(!system.Box(0).is_periodic()) return;

    if(!initialized){
        // Do initial unwrapping
        // Make temp selection from no_jump_ind
        Selection sel(system);
        sel.modify(no_jump_ind);

        // Do unwrapping if more than 1 atom and distance >0
        if(no_jump_ind.size()>1 && unwrap_d>=0){
            LOG()->info("Initial unwrapping of atoms with jump removal...");
            if(unwrap_d==0){
                // Auto find distance
                unwrap_d = 0.2;

                // Find minimal box extent in needed dimensions
                float min_extent = 1e20;
                for(int i=0;i<3;++i)
                    if(dims(i))
                        if(sel.Box().extent(i)<min_extent)
                            min_extent = sel.Box().extent(i);

                while(sel.unwrap_bonds(unwrap_d,leading_index,dims)>1){
                    LOG()->info("Cutoff {} too small for unwrapping.", unwrap_d);
                    unwrap_d *= 2.0;
                    LOG()->info("Trying {}...", unwrap_d);
                    if(unwrap_d > 0.5*min_extent){
                        LOG()->warn("Reached cutoff > 0.5 of box extents!\n"
                                "Selection is likely to consist of disconnected parts.\n"
                                "Continuing as is.");
                        break;
                    }
                }
            } else {
                // Unwrap with given distance
                sel.unwrap_bonds(unwrap_d,leading_index,dims);
            }
            LOG()->info("Unwrapping done.");
        }

        // Save reference coordinates
        no_jump_ref.resize(3,sel.size());
        for(int i=0;i<sel.size();++i){
            no_jump_ref.col(i) = sel.XYZ(i,0);
        }                

        LOG()->info("Will remove jumps for {} atoms", sel.size());

        initialized = true;

    } else { // For other frames, not first

        int ind;
        for(int i=0;i<no_jump_ind.size();++i){
            ind = no_jump_ind[i];
            // Get image closest to running reference
            system.XYZ(ind,0) = system.Box(0).get_closest_image(system.XYZ(ind,0),
                                                                no_jump_ref.col(i),
                                                                dims);
            // Update running reference
            no_jump_ref.col(i) = system.XYZ(ind,0);
        }

    }
}
