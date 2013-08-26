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
