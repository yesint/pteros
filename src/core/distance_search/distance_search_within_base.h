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

#ifndef DISTANCE_SEARCH_WITHIN_BASE_H_INCLUDED
#define DISTANCE_SEARCH_WITHIN_BASE_H_INCLUDED

#include "distance_search_base.h"
#include "atomic_wrapper.h"

namespace pteros {       

class Distance_search_within_base: public Distance_search_base {
protected:
    // Array of atomic bools for used source points
    std::vector<atomic_wrapper<bool>> used;
    void do_search(int sel_size);
    void do_part(int dim, int _b, int _e);
    // Pointer to source selection
    Selection* src_ptr;
    void search_in_pair_of_cells(int sx, int sy, int sz, int tx, int ty, int tz, bool is_periodic);
    void used_to_result(std::vector<int>& res, bool include_self,
                        const Selection &src, const Selection &target);
};

}

#endif
