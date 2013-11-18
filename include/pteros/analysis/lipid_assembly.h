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

#include <fstream>

namespace pteros {

struct Local_curvature {
    // Coordinates of point
    Eigen::Vector3f point;
    Eigen::Vector3f surface_normal;
    Eigen::Vector2f principal_curvatures;
    Eigen::Matrix<float,3,2> principal_directions;
    float gaussian_curvature;
    float mean_curvature;
    float fit_rms;
    float roughness;

    Local_curvature():
        point(Eigen::Vector3f::Zero()),
        surface_normal(Eigen::Vector3f::Zero()),
        principal_curvatures(Eigen::Vector2f::Zero()),
        principal_directions(Eigen::Matrix<float,3,2>::Zero()),
        gaussian_curvature(0.0),
        mean_curvature(0.0),
        fit_rms(0.0),
        roughness(0.0) {}
};

class Lipid_assembly {
public:
    Lipid_assembly(){};
    Lipid_assembly(System& system,
                   std::string lipids_sel,
                   std::string heads_sel,
                   std::string tails_sel,
                   float d = 2.0,
                   int Nsmooth = 1);
    void create(System& system,
                std::string lipids_sel,
                std::string heads_sel,
                std::string tails_sel,
                float d = 2.0,
                int Nsmooth = 1);

    void compute(int frame);
    void write_output();
    Local_curvature weighted_curvature_in_point(Eigen::Vector3f &point);

    int num_lipids(){ return heads.size(); }
    Local_curvature head_curvature(int i){ return head_props[i]; }
    Local_curvature tail_curvature(int i){ return tail_props[i]; }

protected:
    float dist;
    int Nsm;
    // Arrays of lipid properties
    std::vector<Local_curvature> head_props;
    std::vector<Local_curvature> tail_props;
    Selection head_markers;
    Selection tail_markers;
    std::vector<Selection> heads;
    std::vector<Selection> tails;
    Selection lipids;
    int extra_frames; // Counter of created extra aux frames

    void create_markers(std::vector<Selection>& lipids, Selection& markers);
    void get_local_curvature(Selection& surf_spot, Selection* tail_spot,
                             Local_curvature& prop,
                             bool force_flip = false);
    void compute_surface(Selection& markers,
                         std::vector<Selection>& tails,
                         float d, int Nsm,
                         std::vector<Local_curvature>& props,
                         bool force_flip = false);

    void write_vmd_arrows(Selection& markers,
                                          const std::vector<Local_curvature>& props,
                                          std::string fname);

    void write_mean_curvature_as_pdb(Selection& markers,
                                                     const std::vector<Local_curvature>& props,
                                                     std::string fname);
};


}

#endif
