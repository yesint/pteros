/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2023, Semen Yesylevskyy
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


#include "pteros/analysis/jump_remover.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;


JumpRemover::JumpRemover():
    sys_ptr(nullptr),
    dims(fullPBC),
    unwrap_d(0),
    pbc_atom(0),
    initialized(false)
{ }

void JumpRemover::add_atoms(const Selection &sel)
{
    if(!sys_ptr){
        // This is the first invocation
        sys_ptr = sel.get_system();
    } else if(sel.get_system() != sys_ptr){
        throw PterosError("Selections for JumpRemover have to be from the same system!");
    }

    copy(sel.index_begin(),sel.index_end(),back_inserter(no_jump_ind));

    // Only keep unique atoms!
    sort(no_jump_ind.begin(),no_jump_ind.end());
    auto it = unique(no_jump_ind.begin(), no_jump_ind.end());
    no_jump_ind.resize(it - no_jump_ind.begin());
}

void JumpRemover::set_pbc(Array3i_const_ref pbc)
{
    dims = pbc;
    if(dims.sum()==0) LOG()->warn("No periodic dimensions, jump removing disabled.");
}

void JumpRemover::set_unwrap_dist(float d)
{
    unwrap_d = d;
}

void JumpRemover::set_pbc_atom(int ind)
{
    pbc_atom = ind;
}

void JumpRemover::remove_jumps(int fr){
    // Exit immediately if no atoms or no valid dimensions
    // If not periodic also do nothing
    if(no_jump_ind.empty() || dims.sum()==0 || !sys_ptr->box(fr).is_periodic()) return;

    if(!initialized){
        // Do initial unwrapping
        // Make temp selection from no_jump_ind
        Selection sel(*sys_ptr,no_jump_ind);
        sel.set_frame(fr);

        // Do unwrapping if more than 1 atom and distance >=0
        if(no_jump_ind.size()>1 && unwrap_d>=0){
            if(unwrap_d==0){
                if(sel.get_system()->force_field_ready()){
                    // Use topology
                    LOG()->info("Unwrapping using provided topology...");
                    sel.unwrap_bonds(0,dims,pbc_atom);
                } else {
                    // Auto find distance
                    unwrap_d = 0.2;
                    LOG()->info("Trying unwrapping for jump remover, cutoff {}...",unwrap_d);

                    // Find minimal box extent in needed dimensions
                    float min_extent = 1e20;
                    for(int i=0;i<3;++i)
                        if(dims(i))
                            if(sel.box().extent(i)<min_extent)
                                min_extent = sel.box().extent(i);

                    while(sel.unwrap_bonds(unwrap_d,dims,pbc_atom)>1){
                        LOG()->info("Cutoff {} is too small, trying {}...", unwrap_d, 2.0*unwrap_d);
                        unwrap_d *= 2.0;
                        if(unwrap_d > 8.0){
                            LOG()->warn("Cutoff becomes too large! Beware huge memory usage!");
                        }
                        if(unwrap_d > 0.5*min_extent){
                            LOG()->warn("Reached cutoff {} > 0.5 of box extents!\n"
                                        "Selection is likely to consist of disconnected parts.\n"
                                        "Continuing as is.",unwrap_d);
                            break;
                        }
                    }
                    LOG()->info("Unwrapping done at cutoff {}",unwrap_d);
                }
            } else {
                // Unwrap with given distance
                LOG()->info("Unwrapping for jump remover, fixed cutoff {}",unwrap_d);
                sel.unwrap_bonds(unwrap_d,dims,pbc_atom);
            }
        }

        // Save reference coordinates
        no_jump_ref = sel.get_xyz(fr);

        LOG()->info("Will remove jumps for {} atoms", sel.size());

        initialized = true;

    } else { // For other frames, not first

        auto const& box = sys_ptr->box(fr);
        for(size_t i=0;i<no_jump_ind.size();++i){
            int ind = no_jump_ind[i];
            // Get image closest to running reference
            sys_ptr->xyz(ind,fr) = box.closest_image(sys_ptr->xyz(ind,fr),
                                                                  no_jump_ref.col(i),
                                                                  dims);
            // Update running reference
            no_jump_ref.col(i) = sys_ptr->xyz(ind,fr);
        }

    }
}
