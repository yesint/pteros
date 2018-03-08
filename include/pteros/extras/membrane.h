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

    Eigen::Vector3f get_mid_xyz() const {return mid_sel.xyz(0);}
    Eigen::Vector3f get_head_xyz() const {return head_sel.xyz(0);}
    Eigen::Vector3f get_tail_xyz() const {return tail_sel.xyz(0);}

    Selection whole_sel;
    Selection head_sel;
    Selection tail_sel;
    Selection mid_sel;
    // Indexes of carbons in tails for order parameter
    std::vector<std::vector<int>> tail_carbon_indexes;

    std::string name;
    Eigen::Vector3f normal;
    Eigen::Vector3f smoothed_mid_xyz;
    float tilt;    
    float quad_fit_rms;
    float area;
    int leaflet;
    float gaussian_curvature;
    float mean_curvature;
    int coord_number;
    std::vector<std::vector<float>> order;

private:    
    // Set current COM coordinates of seletions to their first atoms used as markers
    void set_markers();
    // Restore atomic coords of markers
    void unset_markers();

    Eigen::Vector3f saved_head0, saved_tail0, saved_mid0;

    Selection local_sel;
};


struct Splay_pair {
    int lip1;
    int lip2;
    float splay;
};

class Membrane {
public:
    Membrane(System *sys, const std::vector<Lipid_descr>& species);
    void compute_properties(float d,
                            Vector3f_const_ref external_normal = Eigen::Vector3f::Zero());
    void write_vmd_arrows(const std::string& fname);
    void write_smoothed(const std::string &fname);

    int num_lipids(){ return lipids.size(); }
    const Lipid& get_lipid(int i){ return lipids[i]; }

    int num_leaflets(){ return leaflets.size(); }
    const std::vector<int>& get_leaflet(int i){ return leaflets[i]; }

    std::vector<Lipid> lipids;
    std::vector<Splay_pair> splay;
    std::vector<std::vector<int>> neighbors;
    std::vector<std::vector<int>> leaflets;

private:
    System* system;
    std::vector<Lipid_descr> lipid_species;
    std::vector<Selection> leaflets_sel;
    std::shared_ptr<spdlog::logger> log;
    std::unordered_map<int,int> index_map;
    // Dynamic properties
    std::vector<Eigen::Vector2i> neighbor_pairs;
};

}

