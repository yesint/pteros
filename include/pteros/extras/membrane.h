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
};



class Lipid {
public:
    Lipid(const Selection& sel, const Lipid_descr& descr);

    void set_leaflet(int i){leaflet = i;}
    // Set current COM coordinates of seletions to their first atoms used as markers
    void set_markers();
    // Restore atomic coords of markers
    void unset_markers();
    Selection get_mid_marker() {return mid_sel(0,0);}
    Selection get_head_marker() {return head_sel(0,0);}
    Selection get_tail_marker() {return tail_sel(0,0);}

    Eigen::Vector3f get_mid_xyz() {return mid_sel.XYZ(0);}
    Eigen::Vector3f get_head_xyz() {return head_sel.XYZ(0);}
    Eigen::Vector3f get_tail_xyz() {return tail_sel.XYZ(0);}

    const Selection& get_mid_sel() {return mid_sel;}

    Eigen::Vector3f normal;
    float angle;
private:
    std::string name;
    Selection whole;
    Selection head_sel;
    Selection tail_sel;
    Selection mid_sel;
    int leaflet;
    Eigen::Vector3f saved_head0, saved_tail0, saved_mid0;
};


class Membrane {
public:
    Membrane(System *sys, const std::vector<Lipid_descr>& species);
    void compute_properties(float d, Vector3f_const_ref external_normal = Eigen::Vector3f::Zero());
    void write_vmd_arrows(const std::string& fname);
private:
    System* system;
    std::vector<Lipid_descr> lipid_species;
    std::vector<Lipid> lipids;
    std::vector<std::vector<int>> leaflets;
    std::vector<Selection> leaflets_sel;

    std::shared_ptr<spdlog::logger> log;

    std::unordered_map<int,int> resindex_map;
};

}
