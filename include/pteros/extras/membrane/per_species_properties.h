#pragma once

#include <Eigen/Core>
#include "pteros/core/utilities.h"
#include "lipid_molecule.h"

namespace pteros {

class LipidMembrane;

class PerSpeciesProperties {
public:
    PerSpeciesProperties(const LipidSpecies *s_ptr, const LipidMembrane *m_ptr);

    double count; // number of lipids of this species. float to avoid overflow.
    // Area
    Histogram area_hist;
    MeanStdAccumulator area;
    // Tilt
    Histogram tilt_hist;
    MeanStdAccumulator tilt;
    // Coordination number
    MeanStdAccumulator coord_number;
    // Trans dihedrals ratio
    MeanStdAccumulator trans_dihedrals_ratio;
    // Curvature
    MeanStdAccumulator gaussian_curvature;
    MeanStdAccumulator mean_curvature;
    Histogram mean_curv_hist;
    Histogram gauss_curv_hist;

    // Order parameters for each tails
    std::vector<Eigen::ArrayXf> order;

    // Histogram of order parameter in normal direction
    // Zero is position of the surf marker for each lipid,
    // Positive values go towards water phase,
    // negative - towards membrane center
    //Histogram order_hist;

    // Abundance of neighboring species
    std::map<std::string,float> around;

    // Called at each lipid on each frame
    void add_data(const LipidMolecule& lip);

    // Called at the end
    void post_process(double num_frames);

    // Returns summary as a string
    std::string summary();

    // Save order to file
    void save_order_to_file(const std::string& fname);
    void save_around_to_file(const std::string& fname);

    int num_tails;
private:    
    const LipidMembrane* membr_ptr;
    const LipidSpecies* sp_ptr;
};

}
