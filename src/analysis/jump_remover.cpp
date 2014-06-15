/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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

using namespace std;
using namespace pteros;


void Jump_remover::add_no_jump_atoms(const Selection &sel)
{    
    int ind;
    int n = sel.get_system()->num_atoms();
    for(int i=0;i<sel.size();++i){
        ind = sel.Index(i);
        if(ind<0 || ind>=n) throw Pteros_error("Index for jump removal out of range!");
        no_jump_ind.push_back(ind);
    }

    // Only keep unique atoms!
    sort(no_jump_ind.begin(),no_jump_ind.end());
    auto it = unique(no_jump_ind.begin(), no_jump_ind.end());
    no_jump_ind.resize( it - no_jump_ind.begin() );
}

void Jump_remover::remove_jumps(System& system, const Frame_info &info){
    // Exit immediately if no atoms
    if(no_jump_ind.empty()) return;        
    // If not periodic also do nothing
    if(!system.Box(0).is_periodic()) return;

    if(info.valid_frame==0){
        // Do initial unwrapping
        // Make temp selection from no_jump_ind
        Selection sel(system);
        sel.modify(no_jump_ind);

        // Do unwrapping if more than 1 atom
        if(no_jump_ind.size()>1){
            cout << "Initial unwrapping of atoms with jump removal..." << endl;
            float cutoff = 0.2;
            float min_extent = system.Box(0).extents().minCoeff();
            while(true){
                try{
                    sel.unwrap_bonds(cutoff);
                }catch(Pteros_error){
                    cout << "Cutoff " << cutoff << " too small for unwrapping. ";
                    cutoff *= 2.0;
                    cout << "Trying " << cutoff << "..." <<endl;
                    if(cutoff > 0.5*min_extent){
                        cout << "Reached cutoff > 0.5 of box extents!\n"
                                "Selection is likely to consist of disconnected parts.\n"
                                "Continuing as is." << endl;
                        break;
                    }
                    continue;
                }
                // If we are here unwrapping is successfull
                break;
            }
            cout << "Unwrapping done." << endl;
        }

        // Save reference coordinates
        no_jump_ref.resize(3,sel.size());
        for(int i=0;i<sel.size();++i){
            no_jump_ref.col(i) = sel.XYZ(i,0);
        }                

        cout << "Removing jumps for " << sel.size() << " atoms" << endl;

    } else { // For other frames, not first

        // remove jumps for every atom of selection. Executed on each step
        int ind;
        for(int i=0;i<no_jump_ind.size();++i){
            ind = no_jump_ind[i];
            // Get image closest to running reference
            system.XYZ(ind,0) = system.Box(0).get_closest_image(system.XYZ(ind,0),
                                                                no_jump_ref.col(i),
                                                                false);
            // Update running reference
            no_jump_ref.col(i) = system.XYZ(ind,0);
        }
    }
}
