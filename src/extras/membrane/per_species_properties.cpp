#include "pteros/extras/membrane/per_species_properties.h"
#include <fstream>
#include "fmt/format.h"
#include "fmt/os.h"
#include "pteros/extras/membrane/lipid_membrane.h"

using namespace pteros;
using namespace std;
using namespace Eigen;

PerSpeciesProperties::PerSpeciesProperties(LipidMembrane *ptr)
{
    membr_ptr = ptr;
    count = 0;

    area_hist.create(0,1.8,100);
    area.fill(0.0);

    tilt_hist.create(0,90,90);
    tilt.fill(0.0);

    coord_number.fill(0.0);

    gaussian_curvature.fill(0.0);
    mean_curvature.fill(0.0);
    mean_curv_hist.create(-0.6,0.6,200);
    gauss_curv_hist.create(-0.3,0.3,200);

    trans_dihedrals_ratio.fill(0.0);
    order_initialized = false;

    // Initialize around
    for(const auto& sp: membr_ptr->species){
        around[sp.name] = 0;
    }

    num_tails = 0;
}


void accumulate_statistics(float val, Eigen::Vector2f& storage){
    storage[0] += val; // mean
    storage[1] += val*val; // std
}


void PerSpeciesProperties::add_data(const LipidMolecule &lip){
    ++count;
    // Area
    area_hist.add(lip.area);
    accumulate_statistics(lip.area,area);

    // Tilt
    tilt_hist.add(lip.tilt);
    accumulate_statistics(lip.tilt,tilt);

    // Coordination number
    accumulate_statistics(lip.coord_number, coord_number);

    // Curvatures
    accumulate_statistics(lip.mean_curvature, mean_curvature);
    accumulate_statistics(lip.gaussian_curvature, gaussian_curvature);
    mean_curv_hist.add(lip.mean_curvature);
    gauss_curv_hist.add(lip.gaussian_curvature);

    // Tail stats
    if(!order_initialized && lip.tails.size()){
        num_tails = lip.tails.size();

        // This is very first invocation, so resize order arrays properly
        order.resize(lip.tails.size());
        for(int i=0;i<order.size();++i){
            order[i].resize(lip.tails[i].order.size());
            order[i].fill(0.0); // Init to zeros
        }

        order_initialized = true;
    } // Initialization


    // Order
    for(int i=0;i<lip.tails.size();++i){
        order[i] += lip.tails[i].order;
    }
    // Trans dihedrals
    for(const auto& t: lip.tails){
        float ratio = (t.dihedrals > M_PI_2).count()/float(t.dihedrals.size());
        accumulate_statistics(ratio,trans_dihedrals_ratio);
    }

    // Lipid surrounding
    for(int i: lip.neib){
        //cout << i << " " << lip.membr_ptr->lipids.size() << endl;
        around[lip.membr_ptr->lipids[i].species_ptr->name] += 1;
    }

}

void PerSpeciesProperties::post_process(float num_frames)
{
    // Skip if no data
    if(count==0 || num_frames==0) return;

    // Compute averages

    // Area
    mean_std_from_accumulated(area,count);
    area_hist.normalize(count);

    // Tilt
    mean_std_from_accumulated(tilt,count);
    tilt_hist.normalize(count);

    // Trans dihedrals
    mean_std_from_accumulated(trans_dihedrals_ratio, count*num_tails); // Note number of tails!

    // Coordination number
    mean_std_from_accumulated(coord_number,count);

    // Curvatures
    mean_std_from_accumulated(mean_curvature,count);
    mean_std_from_accumulated(gaussian_curvature,count);
    mean_curv_hist.normalize(count);
    gauss_curv_hist.normalize(count);

    // Order
    // Average orders for all tails
    for(int i=0;i<order.size();++i) order[i] /= count;

    // At the end set average number of lipids of this kind per frame
    count /= num_frames;

    // Around
    float N = 0;
    for(auto& el: around) N += el.second;
    for(auto& el: around) el.second /= N;
}

string PerSpeciesProperties::summary()
{
    string s;
    if(count>0){
        s += fmt::format("\t\tCount:\t{}\n", count);
        s += fmt::format("\t\tArea:\t{} +/- {} nm2\n", area[0],area[1]);
        s += fmt::format("\t\tTilt:\t{} +/- {} deg\n", rad_to_deg(tilt[0]),rad_to_deg(tilt[1]));
        s += fmt::format("\t\tCoord.N:\t{} +/- {}\n", coord_number[0],coord_number[1]);
        s += fmt::format("\t\tMean.curv.:\t{} +/- {} nm-1\n", mean_curvature[0],mean_curvature[1]);
        s += fmt::format("\t\tGaus.curv.:\t{} +/- {} nm-1\n", gaussian_curvature[0],gaussian_curvature[1]);
        s += fmt::format("\t\tTr.Dih.:\t{} +/- {}\n", trans_dihedrals_ratio[0],trans_dihedrals_ratio[1]);
    } else {
        s += "\t\tNo data\n";
    }
    return s;
}

void PerSpeciesProperties::save_order_to_file(const string &fname)
{
    // Do nothing if no data
    if(count==0 || num_tails==0) return;

    ofstream out(fname);

    // Find the longest tail
    int max_len = 0;
    for(int t=0;t<num_tails;++t)
        if(order[t].size()>max_len) max_len = order[t].size();
    // Header
    fmt::print(out,"#c_num\t");
    for(int t=0;t<num_tails;++t) fmt::print(out,"t{}\t",t);
    fmt::print(out,"\n");
    // Body
    for(int c=0;c<max_len;++c){
        fmt::print(out,"{}\t",c+2);
        for(int t=0;t<num_tails;++t){
            if(c<order[t].size())
                fmt::print(out,"{: .4f}\t",order[t][c]);
            else
                fmt::print(out,"--\t");
        }
        fmt::print(out,"\n");
    }

    out.close();
}

void PerSpeciesProperties::save_around_to_file(const string &fname)
{
    // Do nothing if no data
    if(count==0) return;

    ofstream out(fname);
    for(const auto& sp: membr_ptr->species){
        fmt::print(out,"{}\t{:.4f}\n",sp.name,around[sp.name]);
    }
    out.close();
}
