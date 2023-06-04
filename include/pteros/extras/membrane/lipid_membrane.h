/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2023, Semen Yesylevskyy
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

#include "pteros/extras/membrane/lipid_molecule.h"
#include "pteros/extras/membrane/lipid_group.h"
#include "pteros/extras/membrane/lipid_species.h"
#include "pteros/extras/membrane/lipid_tail.h"

namespace pteros {

void mean_std_from_accumulated(Eigen::Vector2f& storage, float N);

struct InterpolatedPoint {
    InterpolatedPoint():
        normal(Eigen::Vector3f::Zero()),
        mean_curvature(0.0),
        mean_depth(0.0)
    {}

    Eigen::Vector3f normal;
    float mean_curvature;
    float mean_depth;
    std::vector<int> neib_lipids;
    std::vector<float> weights;
};

class LipidMembrane {
public:
    LipidMembrane(const Selection& input_sel,
                  int ngroups,
                  const std::vector<LipidSpecies> &sp_list,
                  const Selection& incl = {},
                  float incl_h_cutoff = 0.5,
                  bool per_carb_normals=false
                  );

    static std::vector<Selection> get_domains(System const& sys, const std::vector<LipidSpecies> &sp_list, float d=0.4);

    void reset_groups();

    void compute_properties(float d = 2.0,
                            float incl_d = 0.5,
                            OrderType order_type = OrderType::SCD_CORR);

    // Interpolating properties in given points.
    void get_interpolation(const Selection& points, std::vector<InterpolatedPoint>& res, Eigen::MatrixXf& normals, float d=3.0);

    /// Returns matrix (n_shells,2)
    /// Each n-th row is averaged curvatures over neigbour shells up to n
    /// col(0) is mean curvature, col(1) is gaussian curvature
    Eigen::MatrixXf get_average_curvatures(int lipid, int n_shells);

    void compute_triangulation();

    void compute_averages();
    void write_averages(const std::string &out_dir=".");
    void write_vmd_visualization(const std::string& path=".");

    std::vector<LipidMolecule> lipids;
    std::vector<LipidGroup> groups;
    std::vector<LipidSpecies> species;

    bool per_carbon_normals;

private:
    const Selection* input_sel_ptr;
    std::shared_ptr<spdlog::logger> log;

    Selection all_surf_sel;
    System surf_sys;

    Selection inclusion;
    float inclusion_h_cutoff;

    Selection all_tails_sel;

    void register_lipid_species(const LipidSpecies &sp);
    void get_initial_normals();
    void create_lipid_patches(const std::vector<Eigen::Vector2i> &bon, const std::vector<float> &dist);
    void add_inclusions(float incl_d);
    void update_local_selection(int i);
    void incremental_triangulation();
    void update_transforms(int i, Vector3f_const_ref normal);
    void inclusion_coord_to_surf_coord(int ind);
    void smoothed_markers_to_local_coord(int ind);
};


} // namespace pteros




