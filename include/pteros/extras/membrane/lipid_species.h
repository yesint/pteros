#pragma once

#include "lipid_tail_descr.h"

namespace pteros {

// Description of one particular lipid species
class LipidSpecies {
    friend class LipidMembrane;
    friend class LipidMolecule;
public:
    LipidSpecies(const std::string& _name,
                 const std::string& _whole_sel_str,
                 const std::string& _head_sel_str,
                 const std::string& _surf_marker_str,
                 const std::vector<std::string>& _tails_descr_strings
                 );
    // Symbolic name of lipid ("POPC")
    std::string name;
    // Selection for the whole lipids ("resname POPC")
    std::string whole_sel_str;
    // Selection for the head marker ("name P")
    std::string head_sel_str;
    // Selection for the surface marker ("name P")
    std::string surf_marker_str;
    // Description strigs for the tails
    std::vector<std::string> tails_descr_strings;
private:
    // Selection for the tail marker. Deduced from tails
    std::string tail_sel_str;
    // Lipid tails
    std::vector<LipidTailDescr> tails_descr;
    // Initialize
    void init(const Selection &lip_sel);
};

}
