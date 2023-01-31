#include "pteros/extras/membrane/lipid_group.h"
#include "pteros/extras/membrane/lipid_membrane.h"
#include "fmt/os.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

LipidGroup::LipidGroup(LipidMembrane *ptr, int id){
    gr_id = id;
    membr_ptr = ptr;
    num_lipids = 0;
    num_frames = 0;
    trans_dihedrals_ratio.reset();
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
    for(auto& it: species_properties){
        num_lipids += it.second.count;        
        trans_dihedrals_ratio.append(it.second.trans_dihedrals_ratio);
    }

    num_lipids = (num_frames) ? num_lipids/num_frames : 0;

    // Compute averages per lipid species
    for(auto& it: species_properties){
        it.second.post_process(num_frames);
    }
}

string LipidGroup::summary() const
{
    string s;
    s += fmt::format("Group #{}:\n",gr_id);
    s += fmt::format("\tNum.lip.:\t{}\n",num_lipids);

    if(num_lipids>0){
        auto [mean,std] = trans_dihedrals_ratio.get_mean_std();
        s += fmt::format("\tTrans.Dih.:\t{:>8.3g} ± {:<8.3g}\n", mean,std);
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

string LipidGroup::properties_table() const
{
    string s;
    s += "Species\tabund%\tTrDih\tTrDihErr\n";
    for(auto& sp: membr_ptr->species){
        s += fmt::format("{}", sp.name);
        const auto& prop = species_properties.at(sp.name);
        s += fmt::format("\t{: .4f}", 100.0*prop.count/float(num_lipids));
        auto [mean,std] = prop.trans_dihedrals_ratio.get_mean_std();
        s += fmt::format("\t{: .4f}\t{: .4f}", mean,std);
        s += "\n";
    }
    return s;
}

void LipidGroup::save_properties_table(const std::filesystem::path &out_dir) const
{
    auto out = fmt::output_file((out_dir / fmt::format("gr{}_properties.dat",get_id())).native());
    out.print(properties_table());
    out.close();
}

void LipidGroup::save_per_species_properties(const filesystem::path &out_dir) const
{
    for(auto& [sp_name,sp]: species_properties){
        if(sp.count>0){
            string file_prefix = out_dir / fmt::format("gr{}_{}_",gr_id,sp_name);
            // Area
            sp.area_hist.save_to_file(file_prefix + "area.dat");
            // Tilt
            sp.tilt_hist.save_to_file(file_prefix + "tilt.dat");
            // Monolayer thickness
            sp.mono_thickness_hist.save_to_file(file_prefix + "mono_thickness.dat");
            // Curvature
            sp.mean_curv_hist.save_to_file(file_prefix + "mean_curv.dat");
            sp.gauss_curv_hist.save_to_file(file_prefix + "gauss_curv.dat");
            // Order
            sp.save_order_to_file(file_prefix + "order.dat");
            // Output order histogram
            //sp.second.order_hist.save_to_file(file_prefix+"order_hist.dat");
            // Around
            sp.save_around_to_file(file_prefix + "around.dat");
        }
    }
}
