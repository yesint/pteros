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

#ifndef DISTANCE_SEARCH_CONTACTS_2SEL_H_INCLUDED
#define DISTANCE_SEARCH_CONTACTS_2SEL_H_INCLUDED

#include "distance_search_contacts.h"

namespace pteros {       

class Distance_search_contacts_2sel: public Distance_search_contacts {
public:
    Distance_search_contacts_2sel(float d, const Selection& sel1, const Selection& sel2,
                                  std::vector<Eigen::Vector2i>& pairs,
                                  bool absolute_index = false,
                                  bool periodic = false,
                                  std::vector<float>* dist_vec = nullptr);
protected:
    void create_grids(const Selection &sel1, const Selection &sel2);

    void do_part(int dim, int _b, int _e,
                 std::deque<Eigen::Vector2i>& bon,
                 std::deque<float>* dist_vec) override;
};

}

#endif
