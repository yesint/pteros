#pragma once

#include "lipid_tail_descr.h"

namespace pteros {

// Description of one particular lipid species
class LipidSpecies {
    friend class LipidMembrane;
    friend class LipidMolecule;
    friend class PerSpeciesProperties;
public:
    LipidSpecies(const std::string& _name,
                 const std::string& _whole_sel_str,
                 const std::string& _head_subsel_str,
                 const std::string& _tail_subsel_str,
                 const std::string& _surf_subsel_str,
                 const std::vector<std::string>& _tails_descr_strings
                 );
    // Symbolic name of lipid ("POPC")
    std::string name;
    // Selection for the whole lipids ("resname POPC")
    std::string whole_sel_str;
    // Selection for the head marker ("name P")
    std::string head_subsel_str;
    // Selection for the tail marker (last atoms of tails typically).
    std::string tail_subsel_str;
    // Selection for the surface marker ("name P")
    std::string surf_subsel_str;
    // Description strigs for the tails
    std::vector<std::string> tails_descr_strings;
private:
    // Lipid tails
    std::vector<LipidTailDescr> tails_descr;
    // Initialize
    void init(const Selection &lip_sel);
};

}
