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


#ifndef LIPID_ASSEMBLY_H
#define LIPID_ASSEMBLY_H

#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/grid_search.h"

namespace pteros {

class Lipid_assembly {
public:
    Lipid_assembly(){};
    Lipid_assembly(Selection& sel, std::string head_marker_atom, float d = 2.0);
    void create(Selection& sel, std::string head_marker_atom, float d = 2.0, float bilayer_cutoff = 0.15);
protected:    
    Selection* source_ptr;
    std::vector<Eigen::Vector3f> surface_normals;
    std::vector<float> surface_mean_angle;
    std::vector<float> surface_curvature;
};

/*
float point_in_membrane(Eigen::Vector3f& point, Selection& head_markers, float d,
                        const Eigen::Vector3i& pbc_dims = Eigen::Vector3i::Ones());
*/
}

#endif
