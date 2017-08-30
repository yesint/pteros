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

#ifndef DISTANCE_SEARCH_WITHIN_INCLUDED
#define DISTANCE_SEARCH_WITHIN_INCLUDED

#include "pteros/core/selection.h"

namespace pteros {       

/// Class for within searching
class Distance_search_within {
public:
    Distance_search_within();

    Distance_search_within(float d,
                           const Selection& src,
                           bool absolute_index = false,
                           bool periodic = false);

    virtual ~Distance_search_within();

    void setup(float d,
               const Selection& src,
               bool absolute_index = false,
               bool periodic = false);

    void search_within(Vector3f_const_ref coord,
                       std::vector<int> &res);

    void search_within(const Selection& target,
                       std::vector<int> &res,
                       bool include_self=true);

private:
    class Distance_search_within_impl;
    std::unique_ptr<Distance_search_within_impl> p;
};

}

#endif
