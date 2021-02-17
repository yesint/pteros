/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2021, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/


#pragma once

#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include "pteros/core/selection.h"

namespace pteros {

/// Convert selection to OpenBabel molecule
void selection_to_obmol(const Selection& sel, OpenBabel::OBMol &mol, bool babel_bonds = true);

}
