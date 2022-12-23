#include "pteros/extras/membrane/lipid_species.h"

using namespace pteros;

LipidSpecies::LipidSpecies(const std::string &_name,
                           const std::string &_whole_sel_str,
                           const std::string &_head_sel_str,
                           const std::string &_tail_sel_str,
                           const std::string &_surf_marker_str,
                           const std::vector<std::string> &_tails_descr_strings):
    name(_name),
    whole_sel_str(_whole_sel_str),
    head_subsel_str(_head_sel_str),
    tail_subsel_str(_tail_sel_str),
    surf_subsel_str(_surf_marker_str),
    tails_descr_strings(_tails_descr_strings)
{}

void LipidSpecies::init(const Selection& lip_sel)
{
    tails_descr.resize(tails_descr_strings.size());
    for(size_t i=0; i<tails_descr.size(); ++i){
        tails_descr[i].init(lip_sel,tails_descr_strings[i]);
    }
}
