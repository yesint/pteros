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

class Local_properties {
public:
    // Coordinates of smoothed point
    Eigen::Vector3f smoothed;
    Eigen::Vector3f surface_normal;
    Eigen::Vector2f principal_curvatures;
    Eigen::Matrix<float,3,2> principal_directions;
    float gaussian_curvature;
    float mean_curvature;
    float fit_rms;
};

class Lipid_assembly {
public:
    Lipid_assembly(){};
    Lipid_assembly(System& system,
                   std::string lipids_sel,
                   std::string heads_sel,
                   std::string tails_sel,
                   float d = 2.0,
                   int Nsm = 1);
    void create(System& system,
                std::string lipids_sel,
                std::string heads_sel,
                std::string tails_sel,
                float d = 2.0,
                int Nsm = 1);

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
