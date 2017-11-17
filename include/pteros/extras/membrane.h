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
    std::string whole;
    std::string head_marker;
    std::string tail_marker;
    std::string position_marker;
};



class Lipid {
public:
    Lipid(const Selection& sel, const Lipid_descr& descr);
    const Selection& get_position_sel() {return position_marker;}
    void set_leaflet(int i){leaflet = i;}
private:
    std::string name;
    Selection whole;
    Selection head_marker;
    Selection tail_marker;
    Selection position_marker;
    int leaflet;
};



class Membrane {
public:
    Membrane(System *sys, const std::vector<Lipid_descr>& species);
private:
    System* system;
    std::vector<Lipid_descr> lipid_species;
    std::vector<Lipid> lipids;
    std::vector<Selection> leaflets;

    std::shared_ptr<spdlog::logger> log;

    std::unordered_map<int,int> resindex_map;
};

}
