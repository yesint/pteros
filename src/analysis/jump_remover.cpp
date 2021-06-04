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


#include "pteros/analysis/jump_remover.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;


JumpRemover::JumpRemover():
    dims(fullPBC),
    unwrap_d(0),
    pbc_atom(0),
    initialized(false)
{ }

void JumpRemover::add_atoms(const Selection &sel)
{    
    int ind;
    int n = sel.get_system()->num_atoms();
    for(int i=0;i<sel.size();++i){
        ind = sel.index(i);
        if(ind<0 || ind>=n) throw PterosError("Index {} for jump removal out of range (0:{})!",ind,n-1);
        no_jump_ind.push_back(ind);
    }

    // Only keep unique atoms!
    sort(no_jump_ind.begin(),no_jump_ind.end());
    auto it = unique(no_jump_ind.begin(), no_jump_ind.end());
    no_jump_ind.resize( it - no_jump_ind.begin() );
}

void JumpRemover::set_pbc(Array3i_const_ref pbc)
{
    dims = pbc;
    if(dims.sum()==0) LOG()->warn("No periodic dimensions, skipping jump removing.");
}

void JumpRemover::set_unwrap_dist(float d)
{
    unwrap_d = d;
}

void JumpRemover::set_pbc_atom(int ind)
{
    pbc_atom = ind;
}

void JumpRemover::remove_jumps(System& system){
    // Exit immediately if no atoms or no valid dimensions
    // If not periodic also do nothing
    if(no_jump_ind.empty() || dims.sum()==0 || !system.box(0).is_periodic()) return;

    if(!initialized){
        // Do initial unwrapping
        // Make temp selection from no_jump_ind
        Selection sel(system,no_jump_ind);

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
        no_jump_ref.resize(3,sel.size());
        for(int i=0;i<sel.size();++i){
            no_jump_ref.col(i) = sel.xyz(i,0);
        }                

        LOG()->info("Will remove jumps for {} atoms", sel.size());

        initialized = true;

    } else { // For other frames, not first

        int ind;
        for(int i=0;i<no_jump_ind.size();++i){
            ind = no_jump_ind[i];
            // Get image closest to running reference
            system.xyz(ind,0) = system.box(0).closest_image(system.xyz(ind,0),
                                                            no_jump_ref.col(i),
                                                            dims);
            // Update running reference
            no_jump_ref.col(i) = system.xyz(ind,0);
        }

    }
}
