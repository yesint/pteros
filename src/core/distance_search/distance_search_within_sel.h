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

#ifndef DISTANCE_SEARCH_WITHIN_SEL_H_INCLUDED
#define DISTANCE_SEARCH_WITHIN_SEL_H_INCLUDED

#include "distance_search_within_base.h"
#include "atomic_wrapper.h"

namespace pteros {       


class Distance_search_within_sel: public Distance_search_within_base {
public:
    /// Constuctor for very fast immediate search of atoms from src,
    /// which are within given distance from the atoms of target.
    /// Used in internal parsing of within selections.
    /// \warning Returns absolute indexes only!
    Distance_search_within_sel(float d,
                           const Selection& src,
                           const Selection& target,
                           std::vector<int> &res,
                           bool include_self=true,
                           bool periodic = false);
};

}

#endif
