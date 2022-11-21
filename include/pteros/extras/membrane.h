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


struct LipidTail;
class LipidMembrane;

struct LocalPatch {
    int lip_id;
    // Id's and distances correspond to each other
    std::vector<int> neib_id;
    std::vector<float> neib_dist;
    // Local coordinate axes in lab space
    Eigen::Matrix3f axes;
    // Transformations to and from local cordinates
    Eigen::Matrix3f to_lab,to_local;
    // Normal
    Eigen::Vector3f normal;
    // Original central point
    Eigen::Vector3f original_center;
};


// Helper class representing a quadric surface for fitting
// Surface is assumed to be centered around point {0,0,0}
// We fit with polynomial fit = A*x^2 + B*y^2 + C*xy + D*x + E*y + F
class QuadSurface {
public:
    // Compute fitting surface given the points in local basis
    void fit_to_points(const Eigen::MatrixXf& coord);

    // Projects point from the XY plane to the surface
    // existing Z coordinate will be substituted in place by the surface Z value
    inline void project_point_to_surface(Vector3f_ref p){
        p(2) =  evalZ(p(0),p(1));
    }

    // Compute Z point of the surface
    inline float evalZ(float x, float y){
        return    A()*x*x
                + B()*y*y
                + C()*x*y
                + D()*x
                + E()*y
                + F();
    }

    // Computes voronoi tesselation in the local tangent plane
    // Sets vertexes of points in plane making up the area
    // Computes in-plane area
    //
    // TODO: Ability to compute area from multiple atoms for each lipid
    // (sum of areas belonging to current lipid)
    //
    void compute_area();

    void compute_curvature();

    inline float A(){ return quad_coefs[0]; }
    inline float B(){ return quad_coefs[1]; }
    inline float C(){ return quad_coefs[2]; }
    inline float D(){ return quad_coefs[3]; }
    inline float E(){ return quad_coefs[4]; }
    inline float F(){ return quad_coefs[5]; }

//qprivate:
    // Coefficients of the quadric surface A,B,C,D,E,F
    Eigen::Matrix<float,6,1> quad_coefs;
    // Fitted points
    Eigen::MatrixXf fitted_points;
    // Fitted normal (unoriented!)
    Eigen::Vector3f fitted_normal;
    // Vertexes for area calculations
    std::vector<Eigen::Vector3f> area_vertexes;
    // List of neighbour ids
    std::vector<int> neib_id;
    // Computed properties
    float fit_rms;
    float in_plane_area;
    float surf_area;
    float gaussian_curvature;
    float mean_curvature;
};

// Properties of individual lipids
struct LipidProperties {
    Eigen::Vector3f normal;
    float tilt;
    float area;
    int coord_number;
    float gaussian_curvature;
    float mean_curvature;
};


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
    Eigen::Vector3f tail_head_vector;

    // Tails
    std::vector<LipidTail> tails;

    // Instanteneous properties of lipid
    Eigen::Vector3f normal;
    float tilt;
    float gaussian_curvature;
    float mean_curvature;

    // Properties from Voronoi tesselation
    int coord_number;
    float area;
    std::vector<int> neib;

    // Fiting resutls
    Eigen::Vector3f smoothed_mid_xyz;
    float quad_fit_rms;

private:

    int id;
    LipidMembrane* membr_ptr;
    // Set markers to current COM coordinates of marker seletions
    void set_markers();
    void unset_markers();

    // Coordinates of markers
    Eigen::Vector3f head_marker, tail_marker, mid_marker, pos_saved;

    Selection local_sel, local_sel_with_self;

    // Staff related to local patch computations
    LocalPatch patch;
    QuadSurface surf;
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
    Histogram mean_curv_hist;
    Histogram gauss_curv_hist;
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

    void compute_properties(float d = 2.0);

    /// Returns matrix (n_shells,2)
    /// Each n-th row is averaged curvatures over neigbour shells up to n
    /// col(0) is mean curvature, col(1) is gaussian curvature
    Eigen::MatrixXf get_average_curvatures(int lipid, int n_shells);

    void compute_averages();
    void write_averages(std::string path=".");
    void write_vmd_visualization(const std::string& path=".");

    std::vector<LipidMolecule> lipids;
    std::vector<LipidGroup> groups;

    std::vector<std::string> species_names;

private:
    System* system; // Parent system
    std::shared_ptr<spdlog::logger> log;

    Selection all_mid_sel;
};


} // namespace pteros




