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

    // Returns summary as a string
    std::string summary();
    // Return properties table as a string
    std::string properties_table();
    // Saves properties to file
    void save_properties_table_to_file(const std::string& fname);

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
