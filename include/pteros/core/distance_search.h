/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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

#ifndef DISTANCE_SEARCH_INCLUDED
#define DISTANCE_SEARCH_INCLUDED

#include "pteros/core/selection.h"
#include "pteros/core/distance_search_within.h"

namespace pteros {       

void search_contacts(float d,
                     const Selection& sel,
                     std::vector<Eigen::Vector2i> &pairs,
                     bool absolute_index = false,
                     bool periodic = false,
                     std::vector<float> *dist_vec = nullptr);

void search_contacts(float d,
                     const Selection& sel1,
                     const Selection& sel2,
                     std::vector<Eigen::Vector2i>& pairs,
                     bool absolute_index = false,
                     bool periodic = false,
                     std::vector<float>* dist_vec = nullptr);

void search_within(float d,
                   const Selection& src,
                   const Selection& target,
                   std::vector<int> &res,
                   bool include_self=true,
                   bool periodic = false);

}

#endif
