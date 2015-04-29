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

#include "pteros/core/distance_search.h"
#include "distance_search_contacts_1sel.h"
#include "distance_search_contacts_2sel.h"
#include "distance_search_within_sel.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

namespace pteros {

void search_contacts(float d, const Selection &sel, std::vector<Vector2i> &pairs, bool absolute_index, bool periodic, std::vector<float> *dist_vec)
{
    Distance_search_contacts_1sel(d,sel,pairs,absolute_index,periodic,dist_vec);
}


void search_contacts(float d, const Selection &sel1, const Selection &sel2, std::vector<Vector2i> &pairs, bool absolute_index, bool periodic, std::vector<float> *dist_vec)
{
    Distance_search_contacts_2sel(d,sel1,sel2,pairs,absolute_index,periodic,dist_vec);
}


void search_within(float d, const Selection &src, const Selection &target, std::vector<int> &res, bool include_self, bool periodic)
{
    Distance_search_within_sel(d,src,target,res,include_self,periodic);
}

}
