#pragma once

#include <openbabel/mol.h>
#include "pteros/core/selection.h"

namespace pteros {

/// Convert selection to OpenBabel molecule
void selection_to_obmol(const Selection& sel, OpenBabel::OBMol &mol, bool babel_bonds = true);

}
