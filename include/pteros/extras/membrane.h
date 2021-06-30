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
#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>

namespace pteros {


struct LipidSpecies {
    // Symbolic name of lipid
    std::string name;
    // Selection for the whole lipids
    std::string whole_str;
    // Selection for the head marker
    std::string head_marker_str;
    // Selection for the tail marker
    std::string tail_marker_str;
    // Selection for the middle marker
    std::string mid_marker_str;
    // Selections for the carbons of each tail
    std::vector<std::string> tail_carbons_str;
};


class LipidTail;

class LipidMolecule {
public:
    LipidMolecule(const Selection& lip_mol, const LipidSpecies& sp);

    std::string name;
    int id; // Unique ID equal to the index of lipid in the Membrane object
    Selection whole_sel;
    Selection head_marker_sel;
    Selection tail_marker_sel;
    Selection mid_marker_sel;

    // Tails
    std::vector<LipidTail> tails;
    int num_tails() const {return tails.size();}

    // General instanteneous properties
    Eigen::Vector3f normal;
    float tilt_angle;
    float area;
    int coord_number;
    Eigen::Vector3f dipole; // Dipole
    float dipole_proj; // Dipole projected onto the normal

    // Curvature-related instanteneous properties
    Eigen::Vector3f smoothed_mid_xyz;
    float quad_fit_rms;
    float gaussian_curvature;
    float mean_curvature;

private:
    // Set markers to current COM coordinates of marker seletions
    void set_markers();
    void unset_markers();

    // Coordinates of markers
    Eigen::Vector3f head_marker, tail_marker, mid_marker, mid_saved;

    Selection local_sel;
};


class LipidTail {
public:
    LipidTail(const Selection& lipid_sel, const std::string& tail_sel_str);

    void compute(const LipidMolecule& lipid);
    int size() const {return carbon_offsets.size();}

    // Order parameters. Size N-2
    Eigen::VectorXf order;
    // Dihedral angles. Size N-3
    Eigen::VectorXf dihedrals;
private:
    // Relative offsets of carbon atoms indexes in whole lipid selection. Size N.
    Eigen::VectorXi carbon_offsets;
};


class PerSpeciesProperties {
public:
    PerSpeciesProperties();

    Eigen::Vector2i count; // number of lipids (mean,std)
    // Area
    Histogram area_hist;
    Eigen::Vector2f area; // (mean,std)
    // Tilt
    Histogram tilt_hist;
    Eigen::Vector2f tilt; // (mean,std)
    // Trans dihedrals ratio
    Eigen::Vector2f trans_dihedrals_ratio; // (mean,std)
    // Order parameter
    std::vector<Eigen::VectorXf> order; //Sz order parameter identical to "gmx order -szonly"
    bool equal_tails; // If true
};


class LipidGroup {
public:
    // Name of the group
    std::string name;
    // Lipids by ID
    std::vector<int> ids;

    // Per group averages (mean,std)
    Eigen::Vector2i num_lipids;
    Eigen::Vector2f trans_dihedrals_ratio;
};

//===================================================================
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

    // Properties
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
    std::vector<std::vector<float>> tail_dihedrals; //Tail dihedral angles
    Eigen::Vector3f dipole; // Dipole
    float dipole_proj; // Dipole projected onto the normal

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
    Histogram tail_dihedrals;
    std::vector<std::vector<float>> order; //Sz order parameter identical to "gmx order -szonly"
    bool equal_tails;
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
    void write_averages(std::string path="");

    void write_vmd_arrows(const std::string& fname);
    void write_smoothed(const std::string &fname);

    int num_lipids(){ return lipids.size(); }
    const Lipid& get_lipid(int i){ return lipids[i]; }

    // Properties
    std::vector<Lipid> lipids; // All per-lipid properties are inside
    std::vector<Splay_pair> splay;
    std::vector<std::vector<int>> neighbors;    

    // Per group averages
    std::vector<Eigen::Vector3f> trans_dih_ratio;

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




