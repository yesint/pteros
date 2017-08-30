/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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

// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
//
// A DSSP reimplementation

#include "mas.h"
#include "dssp.h"
#include "structure.h"
#include "pteros_dssp_wrapper.h"

int VERBOSE = 0;

namespace pteros {

void dssp_wrapper(pteros::Selection& sel, std::ostream& io){
    MProtein a(sel);
    a.CalculateSecondaryStructure();
    WriteDSSP(a, io);
}

std::string dssp_string(pteros::Selection& sel){
    MProtein a(sel);
    a.CalculateSecondaryStructure();
    // Now form an output string
    return make_DSSP_string(a);
}

}
