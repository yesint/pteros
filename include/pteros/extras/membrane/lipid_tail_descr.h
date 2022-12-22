#pragma once

#include "pteros/core/selection.h"

namespace pteros {

class LipidTailDescr {
public:
    // tail_descr_str has the following format:
    // C1-C2-C3=C4-C5-C6=C7-C8
    // Where "-" is a single bond and "=" is a double bond
    // This is needed for correct order computations
    LipidTailDescr(){}

    void init(const Selection& lipid_sel, const std::string& tail_descr_str);
    int size() const {return c_names.size();}

    std::vector<std::string> c_names;
    std::vector<int> bond_orders;
    // Relative offsets of carbon atoms indexes in whole lipid selection. Size N.
    std::vector<int> c_offsets;
};

}
