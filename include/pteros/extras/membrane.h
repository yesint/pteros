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
class LipidMembrane;

class LipidMolecule {
friend class LipidMembrane;
friend class PerSpeciesProperties;
public:
    LipidMolecule(const Selection& lip_mol, const LipidSpecies& sp, int ind, LipidMembrane* parent);

    void add_to_group(int gr);

    std::string name;
    Selection whole_sel;
    Selection head_marker_sel;
    Selection tail_marker_sel;
    Selection mid_marker_sel;

    // Tails
    std::vector<LipidTail> tails;

    // General instanteneous properties
    Eigen::Vector3f normal;
    float tilt;
    float area;
    int coord_number;
    Eigen::Vector3f dipole; // Dipole
    float dipole_proj; // Dipole projected onto the normal

    // Curvature-related instanteneous properties
    Eigen::Vector3f smoothed_mid_xyz;
    float quad_fit_rms;
    float gaussian_curvature;
    float mean_curvature;

    // Indexes of neighbour lipids
    std::vector<int> neib;

private:

    int id;
    LipidMembrane* membr_ptr;
    // Set markers to current COM coordinates of marker seletions
    void set_markers();
    void unset_markers();

    // Coordinates of markers
    Eigen::Vector3f head_marker, tail_marker, mid_marker, pos_saved;

    Selection local_sel;
};


struct LipidTail {
    LipidTail(const Selection& lipid_sel, const std::string& tail_sel_str);
    void compute(const LipidMolecule& lipid);
    int size() const {return carbon_offsets.size();}

    // Order parameters. Size N-2
    Eigen::ArrayXf order;
    // Dihedral angles. Size N-3
    Eigen::ArrayXf dihedrals;
    // Relative offsets of carbon atoms indexes in whole lipid selection. Size N.
    Eigen::VectorXi carbon_offsets;
};


class PerSpeciesProperties {    
public:
    PerSpeciesProperties(LipidMembrane* ptr);

    float count; // number of lipids of this species. float to avoid overflow.
    // Area
    Histogram area_hist;
    Eigen::Vector2f area; // (mean,std)
    // Tilt
    Histogram tilt_hist;
    Eigen::Vector2f tilt; // (mean,std)
    // Total dipole
    Eigen::Vector2f total_dipole; // (mean,std)
    // Projected dipole
    Eigen::Vector2f projected_dipole; // (mean,std)
    // Coordination number
    Eigen::Vector2f coord_number; // (mean,std)
    // Trans dihedrals ratio
    Eigen::Vector2f trans_dihedrals_ratio; // (mean,std)
    // Curvature
    Eigen::Vector2f gaussian_curvature;
    Eigen::Vector2f mean_curvature;
    // Order parameter
    std::vector<Eigen::ArrayXf> order; //Sz order parameter identical to "gmx order -szonly"
    // Abundance of neighboring species
    std::map<std::string,float> around;

    // Called at each lipid on each frame
    void add_data(const LipidMolecule& lip);
    // Called at the end
    void post_process(float num_frames);

    // Returns summary as a string
    std::string summary();

    // Save order to file
    void save_order_to_file(const std::string& fname);
    void save_around_to_file(const std::string& fname);

    int num_tails;
private:
    bool order_initialized;
    LipidMembrane* membr_ptr;
};


class LipidGroup {
public:
    LipidGroup(LipidMembrane* ptr, int id);

    void reset(){ lip_ids.clear(); }
    void add_lipid_id(int i){lip_ids.push_back(i);}
    void process_frame();
    void post_process();

    // Returns summary as a string
    std::string summary();
    std::string properties_table();
    void save_properties_table_to_file(const std::string& fname);

    // Per group averages (mean,std)
    float num_lipids, num_frames;
    Eigen::Vector2f trans_dihedrals_ratio;
    // Per species averages
    std::map<std::string,PerSpeciesProperties> species_properties;

private:
    // Group ID
    int gr_id;
    // Lipids by ID
    std::vector<int> lip_ids;
    // Parent ptr
    LipidMembrane* membr_ptr;
};


class LipidMembrane {
public:
    LipidMembrane(System *sys, const std::vector<LipidSpecies>& species, int ngroups);

    void reset_groups();

    void compute_properties(float d = 2.0,
                            bool use_external_normal = false,
                            Vector3f_const_ref external_pivot = Eigen::Vector3f::Zero(),
                            Vector3i_const_ref external_dist_dim = Eigen::Vector3i::Ones());

    void compute_averages();
    void write_averages(std::string path=".");

    std::vector<LipidMolecule> lipids;
    std::vector<LipidGroup> groups;

    std::vector<std::string> species_names;
private:
    System* system; // Parent system
    std::shared_ptr<spdlog::logger> log;

    Selection all_mid_sel;
};


} // namespace pteros




