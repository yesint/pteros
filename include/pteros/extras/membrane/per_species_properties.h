#pragma once

#include <Eigen/Core>
#include "pteros/core/utilities.h"
#include "lipid_molecule.h"

namespace pteros {

class LipidMembrane;

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
    // Coordination number
    Eigen::Vector2f coord_number; // (mean,std)
    // Trans dihedrals ratio
    Eigen::Vector2f trans_dihedrals_ratio; // (mean,std)
    // Curvature
    Eigen::Vector2f gaussian_curvature;
    Eigen::Vector2f mean_curvature;
    Histogram mean_curv_hist;
    Histogram gauss_curv_hist;

    // Order parameters for each tails
    std::vector<Eigen::ArrayXf> order;

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

}
