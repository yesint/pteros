#include "pteros/extras/membrane/per_species_properties.h"
#include "fmt/format.h"
#include "fmt/os.h"
#include "pteros/extras/membrane/lipid_membrane.h"

using namespace pteros;
using namespace std;
using namespace Eigen;

PerSpeciesProperties::PerSpeciesProperties(const LipidSpecies *s_ptr, const LipidMembrane *m_ptr)
{
    membr_ptr = m_ptr;
    sp_ptr = s_ptr;
    count = 0;

    area_hist.create(0,1.8,100);
    area.reset();

    tilt_hist.create(0,90,90);
    tilt.reset();

    mono_thickness_hist.create(0.0,3.0,100);
    mono_thickness.reset();

    coord_number.reset();

    gaussian_curvature.reset();
    mean_curvature.reset();
    mean_curv_hist.create(-0.6,0.6,200);
    gauss_curv_hist.create(-0.3,0.3,200);

    trans_dihedrals_ratio.reset();

    //order_hist.create(-2.5,0.5,100);
    // Resize order arrays properly

    num_tails = sp_ptr->tails_descr.size();
    order.resize(num_tails);
    for(size_t i=0;i<order.size();++i){
        order[i].resize(sp_ptr->tails_descr[i].size()-2);
        order[i].fill(0.0); // Init to zeros
    }

    // Initialize around
    for(const auto& sp: membr_ptr->species){
        around[sp.name] = 0;
    }
}


void accumulate_statistics(float val, Eigen::Vector2f& storage){
    storage[0] += val; // mean
    storage[1] += val*val; // std
}


void PerSpeciesProperties::add_data(const LipidMolecule &lip){
    ++count;
    // Area
    area_hist.add(lip.area);
    area.add_value(lip.area);

    // Tilt
    tilt_hist.add(lip.tilt);
    tilt.add_value(lip.tilt);

    // Monolayer thickness
    mono_thickness_hist.add(lip.mono_thickness);
    mono_thickness.add_value(lip.mono_thickness);

    // Coordination number
    coord_number.add_value(lip.coord_number);

    // Curvatures    
    mean_curv_hist.add(lip.mean_curvature);
    gauss_curv_hist.add(lip.gaussian_curvature);

    mean_curvature.add_value(lip.mean_curvature);
    gaussian_curvature.add_value(lip.gaussian_curvature);

    // Order
    for(size_t t=0;t<lip.tails.size();++t){
        order[t] += lip.tails[t].order;
        /*
        // Add to order z-hist
        for(long a=0; a<order[t].size(); ++a){
            // Find coordinate of tail atom a and compute vector from this atom to surf marker
            // Project this vector to the normal and compute it length
            int loc_ind = sp_ptr->tails_descr[t].c_offsets[a];
            Vector3f v = lip.whole_sel.xyz(loc_ind) - lip.patch.original_center;
            float zd = v.dot(lip.normal)/lip.normal.dot(lip.normal);
            // Compute approx. volume of the slice taking into account local curvature
            //
            order_hist.add(zd,order[t](a)/vol);
        }
        */
    }    


    // Trans dihedrals
    for(const auto& t: lip.tails){
        float ratio = (t.dihedrals > M_PI_2).count()/float(t.dihedrals.size());
        trans_dihedrals_ratio.add_value(ratio);
    }

    // Lipid surrounding
    for(int i: lip.neib){
        //cout << i << " " << lip.membr_ptr->lipids.size() << endl;
        around[lip.membr_ptr->lipids[i].species_ptr->name] += 1;
    }

}

void PerSpeciesProperties::post_process(double num_frames)
{
    // Skip if no data
    if(count==0 || num_frames==0) return;

    // Compute averages

    // Area
    area_hist.normalize(count);

    // Tilt
    tilt_hist.normalize(count);

    // Thickness
    mono_thickness_hist.normalize(count);

    // Curvatures
    mean_curv_hist.normalize(count);
    gauss_curv_hist.normalize(count);

    // Order
    // Average orders for all tails
    for(size_t i=0;i<order.size();++i) order[i] /= count;
    //order_hist.normalize(count);

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
    Vector2d v;
    if(count>0){
        s += fmt::format("\t\tCount:\t\t{:>8g}\n", count);

        v = area.get_mean_std();
        s += fmt::format("\t\tArea:\t\t{:>8.3g} ± {:<8.3g} nm²\n", v(0),v(1));

        v = tilt.get_mean_std();
        s += fmt::format("\t\tTilt:\t\t{:>8.3g} ± {:<8.3g} deg\n", v(0),v(1));

        v = coord_number.get_mean_std();
        s += fmt::format("\t\tCoord.N:\t{:>8.3g} ± {:<8.3g}\n", v(0),v(1));

        v = mono_thickness.get_mean_std();
        s += fmt::format("\t\tMono.thick.:\t{:>8.3g} ± {:<8.3g} nm\n", v(0),v(1));

        v = mean_curvature.get_mean_std();
        s += fmt::format("\t\tMean.curv.:\t{:>8.3g} ± {:<8.3g} nm⁻¹\n", v(0),v(1));

        v = gaussian_curvature.get_mean_std();
        s += fmt::format("\t\tGaus.curv.:\t{:>8.3g} ± {:<8.3g} nm⁻¹\n", v(0),v(1));

        v = trans_dihedrals_ratio.get_mean_std();
        s += fmt::format("\t\tTrans.Dih.:\t{:>8.3g} ± {:<8.3g}\n", v(0),v(1));
    } else {
        s += "\t\tNo data\n";
    }
    return s;
}

void PerSpeciesProperties::save_order_to_file(const string &fname)
{
    num_tails = order.size();
    // Do nothing if no data
    if(count==0 || num_tails==0) return;

    auto out = fmt::output_file(fname);

    // Find the longest tail
    int max_len = 0;
    for(int t=0;t<num_tails;++t)
        if(order[t].size()>max_len) max_len = order[t].size();
    // Header
    out.print("#c_num\t");
    for(int t=0;t<num_tails;++t) out.print("t{}\t",t);
    out.print("\n");
    // Body
    for(int c=0;c<max_len;++c){
        out.print("{}\t",c+2);
        for(int t=0;t<num_tails;++t){
            if(c<order[t].size())
                out.print("{: .4f}\t",order[t][c]);
            else
                out.print("--\t");
        }
        out.print("\n");
    }

    out.close();    
}

void PerSpeciesProperties::save_around_to_file(const string &fname)
{
    // Do nothing if no data
    if(count==0) return;

    auto out = fmt::output_file(fname);
    for(const auto& sp: membr_ptr->species){
        out.print("{}\t{:.4f}\n",sp.name,around[sp.name]);
    }
    out.close();
}
