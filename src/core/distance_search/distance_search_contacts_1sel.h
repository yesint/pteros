/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/




#pragma once

#include "distance_search_contacts.h"

namespace pteros {       

class DistanceSearchContacts1sel: public DistanceSearchContacts {
public:

    DistanceSearchContacts1sel(float d,
                                  const Selection& sel,
                                  std::vector<Eigen::Vector2i>& res_pairs,
                                  std::vector<float>& res_distances,
                                  bool absolute_index = false,
                                  Vector3i_const_ref pbc = fullPBC);
protected:    

    virtual void search_planned_pair(const PlannedPair& pair,
                                     std::vector<Eigen::Vector2i> &pairs_buffer,
                                     std::vector<float> &distances_buffer) override;
};

}




