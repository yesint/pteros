#pragma once

#include "lipid_tail_descr.h"

namespace pteros {

// Types of order parameter to compute
enum struct OrderType {
    // Sz order parameter identical to gromacs -szonly option
    SZ,
    // Deuterium order parameter computed for ideal H positions for double bonds
    SCD,
    // Deuterium order parameter computed for corrected H positions for double bonds
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3882000/
    SCD_CORR
};


class LipidTail {
public:
    LipidTail(LipidTailDescr* descr);
    // Order parameters. Size N-2
    Eigen::ArrayXf order;

    // Dihedral angles. Size N-3
    Eigen::ArrayXf dihedrals;

    void compute_order_and_dihedrals(const Selection& whole_lipid_sel,
                                     Vector3f_const_ref normal,
                                     OrderType order_type = OrderType::SCD_CORR);
private:
    LipidTailDescr* descr_ptr;
};

}
