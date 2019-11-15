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


#pragma once
#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>

namespace pteros {


struct Lipid_descr {
    std::string name;
    std::string whole_sel_str;
    std::string head_sel_str;
    std::string tail_sel_str;
    std::string mid_sel_str;
    std::vector<std::string> tail_carbon_sels;
};

class Lipid {
    friend class Membrane;
public:
    Lipid(const Selection& sel, const Lipid_descr& descr);

    Selection whole_sel;
    Selection head_sel;
    Selection tail_sel;
    Selection mid_sel;    

    std::string name;
    int group;

    Eigen::Vector3f normal;
    Eigen::Vector3f smoothed_mid_xyz;
    float tilt;    
    float quad_fit_rms;
    float area;
    int leaflet;
    float gaussian_curvature;
    float mean_curvature;
    int coord_number;
    std::vector<std::vector<float>> order; //Sz order parameter identical to "gmx order -szonly"

private:    
    // Set markers to current COM coordinates of seletions
    void set_markers();
    void unset_markers();

    // Coordinates of markers
    Eigen::Vector3f head_marker, tail_marker, mid_marker, mid_saved;

    Selection local_sel;

    // Indexes of carbons in tails for order parameter
    std::vector<std::vector<int>> tail_carbon_indexes;
};


struct Splay_pair {
    int lip1;
    int lip2;
    float splay;
};


struct Average_props_per_type {
    Average_props_per_type();

    int num;
    Histogram area;
    Histogram tilt;
    Histogram gaussian_curvature;
    Histogram mean_curvature;
    Histogram coord_number;
    std::vector<std::vector<float>> order; //Sz order parameter identical to "gmx order -szonly"
};

using Lipid_group = std::map<std::string,Average_props_per_type>;


class Membrane {
public:
    Membrane(System *sys, const std::vector<Lipid_descr>& species,
             int ngroups=1, bool compute_splay=false);


    void compute_properties(float d,
                            bool use_external_normal = false,
                            Vector3f_const_ref external_pivot = Eigen::Vector3f::Zero(),
                            Vector3i_const_ref external_dist_dim = Eigen::Vector3i::Ones());

    void compute_averages();
    void write_averages();

    void write_vmd_arrows(const std::string& fname);
    void write_smoothed(const std::string &fname);

    int num_lipids(){ return lipids.size(); }
    const Lipid& get_lipid(int i){ return lipids[i]; }

    // Properties
    std::vector<Lipid> lipids; // All per-lipid properties are inside
    std::vector<Splay_pair> splay;
    std::vector<std::vector<int>> neighbors;

    // Per group averages


private:
    System* system;
    std::vector<Lipid_descr> lipid_species;    
    std::shared_ptr<spdlog::logger> log;
    std::unordered_map<int,int> index_map;
    // Dynamic properties
    std::vector<Eigen::Vector2i> neighbor_pairs;
    Selection all_mid_sel;
    std::vector<Lipid_group> groups;
    bool do_splay;
};

}

