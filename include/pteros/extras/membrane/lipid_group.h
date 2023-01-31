#pragma once

#include <vector>
#include <string>
#include <map>
#include <Eigen/Core>
#include "per_species_properties.h"

namespace pteros {

class LipidMembrane;

class LipidGroup {
public:
    LipidGroup(LipidMembrane* ptr, int id);

    void reset(){ lip_ids.clear(); }
    void add_lipid_id(int i){lip_ids.push_back(i);}
    void process_frame();
    void post_process();

    int get_id() const {return gr_id;}

    // Returns summary as a string
    std::string summary() const;
    // Return properties table as a string
    std::string properties_table() const;
    // Saves properties to file
    void save_properties_table(const std::filesystem::path &out_dir) const;
    void save_per_species_properties(const std::filesystem::path& out_dir) const;

    // Per group averages (mean,std)
    float num_lipids, num_frames;
    MeanStdAccumulator trans_dihedrals_ratio;
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


}
