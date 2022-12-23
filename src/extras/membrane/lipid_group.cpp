#include "pteros/extras/membrane/lipid_group.h"
#include "pteros/extras/membrane/lipid_membrane.h"
#include <fstream>

using namespace std;
using namespace pteros;
using namespace Eigen;

LipidGroup::LipidGroup(LipidMembrane *ptr, int id){
    gr_id = id;
    membr_ptr = ptr;
    num_lipids = 0;
    num_frames = 0;
    trans_dihedrals_ratio.fill(0.0);
    // Initialize species_properties
    for(const auto& sp: membr_ptr->species){
        species_properties.emplace(sp.name,PerSpeciesProperties(&sp,membr_ptr));
    }
}

void LipidGroup::process_frame()
{
    // Cycle over lipids in this group
    for(int id: lip_ids){
        auto name = membr_ptr->lipids[id].species_ptr->name;
        // Add data to particular scpecies
        species_properties.at(name).add_data(membr_ptr->lipids[id]);
    }

    ++num_frames;
}



void LipidGroup::post_process()
{
    // Collect bulk statistics for the group
    int num_dihedrals = 0.0;
    for(auto& it: species_properties){
        num_lipids += it.second.count;
        num_dihedrals += it.second.count*it.second.num_tails;
        trans_dihedrals_ratio += it.second.trans_dihedrals_ratio;
    }

    // Compute correct averages for bulk properties
    mean_std_from_accumulated(trans_dihedrals_ratio, num_dihedrals);
    num_lipids = (num_frames) ? num_lipids/num_frames : 0;

    // Compute averages per lipid species
    for(auto& it: species_properties){
        it.second.post_process(num_frames);
    }
}

string LipidGroup::summary()
{
    string s;
    s += fmt::format("Group #{}:\n",gr_id);
    s += fmt::format("\tNum.lip.:\t{}\n",num_lipids);

    if(num_lipids>0){
        s += fmt::format("\tTrans.Dih.:\t{:>8.3g} ± {:<8.3g}\n", trans_dihedrals_ratio[0],trans_dihedrals_ratio[1]);
        s += "\tLipid species:\n";
        for(auto& sp: membr_ptr->species){
            s += fmt::format("\t{}:\n", sp.name);
            s += species_properties.at(sp.name).summary();
        }

        // Write table summary
        s += "\n\tProperties table:\n";
        s += properties_table();
    } else {
        s += "\tNo data\n";
    }
    return s;
}

string LipidGroup::properties_table()
{
    string s;
    s += "Species\tabund%\tTrDih\tTrDihErr\n";
    for(auto& sp: membr_ptr->species){
        s += fmt::format("{}", sp.name);
        const auto& prop = species_properties.at(sp.name);
        s += fmt::format("\t{: .4f}", 100.0*prop.count/float(num_lipids));
        s += fmt::format("\t{: .4f}\t{: .4f}", prop.trans_dihedrals_ratio[0],prop.trans_dihedrals_ratio[1]);
        s += "\n";
    }
    return s;
}

void LipidGroup::save_properties_table_to_file(const string &fname)
{
    ofstream out(fname);
    out << properties_table();
    out.close();
}
