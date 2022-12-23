#include "pteros/extras/membrane/lipid_molecule.h"
#include "pteros/extras/membrane/lipid_species.h"
#include "pteros/extras/membrane/lipid_membrane.h"
#include "pteros/core/pteros_error.h"

using namespace pteros;
using namespace Eigen;
using namespace std;

LipidMolecule::LipidMolecule(const Selection& lip_mol, LipidSpecies* sp, int ind, LipidMembrane* parent)
{
    id = ind;
    membr_ptr = parent;
    species_ptr = sp;

    whole_sel = lip_mol;
    head_marker_sel = whole_sel(sp->head_subsel_str);
    tail_marker_sel = whole_sel(sp->tail_subsel_str);
    mid_marker_sel = whole_sel(sp->surf_subsel_str);

    // Initialize tails
    for(int i=0; i<sp->tails_descr.size(); ++i){
        tails.emplace_back(&sp->tails_descr[i]);
    }
}

void LipidMolecule::add_to_group(int gr){
    if(gr<0 || gr>=membr_ptr->groups.size()) throw PterosError("The group should be in the range (0:{}), not {}!",
                                                               membr_ptr->groups.size(),gr);
    membr_ptr->groups[gr].add_lipid_id(id);
}

void LipidMolecule::set_markers()
{
    // Unwrap this lipid with leading index of mid marker
    whole_sel.unwrap(fullPBC, mid_marker_sel.index(0)-whole_sel.index(0));

    // Set markers to COM
    head_marker = head_marker_sel.center(true);
    tail_marker = tail_marker_sel.center(true);
    mid_marker = mid_marker_sel.center(true);

    pos_saved = mid_marker_sel.xyz(0);
    mid_marker_sel.xyz(0) = mid_marker;
    tail_head_vector = head_marker-tail_marker;
}

void LipidMolecule::unset_markers()
{
    mid_marker_sel.xyz(0) = pos_saved;
}
