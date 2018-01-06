/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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

