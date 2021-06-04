#pragma once

#include "pteros/core/selection.h"

namespace pteros {
    void dssp_wrapper(Selection& sel, std::ostream &io);
    std::string dssp_string(pteros::Selection& sel);
}

