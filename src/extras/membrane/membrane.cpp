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

#include "pteros/extras/membrane.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include "voro++.hh"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
    #define M_PI_2 9.86960440108935712052951
#endif

using namespace std;
using namespace pteros;
using namespace Eigen;


string tcl_arrow(Vector3f_const_ref p1, Vector3f_const_ref p2, float r, string color){
    stringstream ss;
    Vector3f p = (p2-p1)*0.8+p1;
    ss << "draw color " << color << endl;

    ss << "draw cylinder \"" << p1.transpose()*10.0 << "\" ";
    ss << "\"" << p.transpose()*10.0 << "\" radius " << r << endl;

    ss << "draw cone \"" << p.transpose()*10.0 << "\" ";
    ss << "\"" << p2.transpose()*10.0 << "\" radius " << r*3.0 << endl;

    return ss.str();
}


//--------------------------------------------


LipidTail::LipidTail(const Selection& lipid_sel, const string &tail_sel_str)
{
    Selection tail_sel = const_cast<Selection&>(lipid_sel)(tail_sel_str);
    int N = tail_sel.size();

    // Set offsets
    carbon_offsets.resize(N);
    for(int i=0;i<N;++i) carbon_offsets[i] = tail_sel.index(i) - lipid_sel.index(0);

    // Initialize order
    order.resize(N-2);
    order.fill(0.0);

    // Initialize dihedrals
    dihedrals.resize(N-3);
    dihedrals.fill(0.0);
}

void LipidTail::compute(const LipidMolecule& lipid)
{
    int N = carbon_offsets.size();
    // Compute order
    for(int at=1; at<N-1; ++at){
        // Vector from at+1 to at-1
        auto coord1 = lipid.whole_sel.xyz(carbon_offsets[at+1]);
        auto coord2 = lipid.whole_sel.xyz(carbon_offsets[at-1]);
        float ang = angle_between_vectors(coord1-coord2,lipid.normal);
        order[at-1] = 1.5*pow(cos(ang),2)-0.5;
    }

    // Compute dihedrals
    for(int at=0; at<N-3; ++at){
        dihedrals[at] = lipid.whole_sel.dihedral(carbon_offsets[at],
                                                 carbon_offsets[at+1],
                                                 carbon_offsets[at+2],
                                                 carbon_offsets[at+3],
                                                 noPBC);
    }
}

LipidMolecule::LipidMolecule(const Selection& lip_mol, const LipidSpecies& sp, int ind, LipidMembrane* parent)
{
    id = ind;
    membr_ptr = parent;
    name = sp.name;
    whole_sel = lip_mol;
    head_marker_sel = whole_sel(sp.head_marker_str);
    tail_marker_sel = whole_sel(sp.tail_marker_str);
    mid_marker_sel = whole_sel(sp.mid_marker_str);
    // Create tails
    tails.reserve(sp.tail_carbons_str.size());
    for(const string& t_str: sp.tail_carbons_str){
        tails.emplace_back(whole_sel, t_str);
    }
    //neib.clear();
}

void LipidMolecule::add_to_group(int gr){
    if(gr<0 || gr>=membr_ptr->groups.size()) throw PterosError("The group should be in the range (0:{}), not {}!",
                                                               membr_ptr->groups.size(),gr);
    membr_ptr->groups[gr].add_lipid_id(id);
}

void LipidMolecule::set_markers()
{
    // Unwrap this lipid with leading index of mid marker
    whole_sel.unwrap(fullPBC, mid_marker_sel.index(0)-whole_sel.index(0));

    // Set markers to COM
    head_marker = head_marker_sel.center();
    tail_marker = tail_marker_sel.center();
    mid_marker = mid_marker_sel.center();

    pos_saved = mid_marker_sel.xyz(0);
    mid_marker_sel.xyz(0) = mid_marker;
    tail_head_vector = head_marker-tail_marker;
}

void LipidMolecule::unset_markers()
{
    mid_marker_sel.xyz(0) = pos_saved;
}

PerSpeciesProperties::PerSpeciesProperties(LipidMembrane *ptr)
{
    membr_ptr = ptr;
    count = 0;

    area_hist.create(0,1.8,100);
    area.fill(0.0);

    tilt_hist.create(0,90,90);
    tilt.fill(0.0);

    total_dipole.fill(0.0);
    projected_dipole.fill(0.0);
    coord_number.fill(0.0);

    gaussian_curvature.fill(0.0);
    mean_curvature.fill(0.0);
    mean_curv_hist.create(-0.6,0.3,100);
    gauss_curv_hist.create(-0.3,0.3,100);

    trans_dihedrals_ratio.fill(0.0);
    order_initialized = false;

    // Initialize around
    for(const auto& sp_name: ptr->species_names){
        around[sp_name] = 0;
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

    // Dipoles
    accumulate_statistics(lip.dipole.norm(), total_dipole);
    accumulate_statistics(lip.dipole_proj, projected_dipole);

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
        // See if the tails are of the same length
        int sz = lip.tails[0].size();
        bool same = true;
        for(int i=1;i<lip.tails.size();++i){
            if(lip.tails[i].size()!=sz){
                same = false;
                break;
            }
        }

        if(!same){
            order.resize(lip.tails.size());
            for(int i=0;i<order.size();++i){
                order[i].resize(lip.tails[i].order.size());
                order[i].fill(0.0); // Init to zeros
            }
        } else {
            // tails are the same add extra slot for average
            order.resize(lip.tails.size()+1);
            for(int i=0;i<order.size();++i){
                order[i].resize(lip.tails[0].order.size()); // Resize also average correctly
                order[i].fill(0.0); // Init to zeros
            }
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
    /*
    for(int i: lip.neib){
        //cout << i << " " << lip.membr_ptr->lipids.size() << endl;
        around[lip.membr_ptr->lipids[i].name] += 1;
    }
    */

}

void mean_std_from_accumulated(Vector2f& storage, float N){
    if(N>0){
        float s1 = storage[0];
        float s2 = storage[1];
        storage[0] = s1/N;
        storage[1] = sqrt(s2/N - s1*s1/N/N);
    } else {
        storage.fill(0.0);
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

    // Dipoles
    mean_std_from_accumulated(total_dipole,count);
    mean_std_from_accumulated(projected_dipole,count);

    // Coordination number
    mean_std_from_accumulated(coord_number,count);

    // Curvatures
    mean_std_from_accumulated(mean_curvature,count);
    mean_std_from_accumulated(gaussian_curvature,count);
    mean_curv_hist.normalize(count);
    gauss_curv_hist.normalize(count);

    // Order
    if(num_tails<order.size()){
        // we have an extra slot - this is for average of the same size tails
        for(int i=0;i<num_tails;++i) order[num_tails] += order[i]/float(num_tails);
    }
    // Average orders for all tails (including the average if present)
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
        s += fmt::format("\t\tDipole:\t{} +/- {} D\n", total_dipole[0],total_dipole[1]);
        s += fmt::format("\t\tDip.proj.:\t{} +/- {} D\n", projected_dipole[0],projected_dipole[1]);
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
    if(num_tails<order.size()){
        // we have an extra slot - this is for average of the same size tails
        // Header
        fmt::print(out,"#c_num\t");
        for(int t=0;t<num_tails;++t) fmt::print(out,"t{}\t",t);
        fmt::print(out,"t_aver\n");

        for(int c=0;c<order[0].size();++c){            
            fmt::print(out,"{}\t",c+2);
            for(int t=0;t<order.size();++t) fmt::print(out,"{: .4f}\t",order[t][c]);
            fmt::print(out,"\n");
        }
    } else {
        // Tails are of different length
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
    }
    out.close();
}

void PerSpeciesProperties::save_around_to_file(const string &fname)
{
    // Do nothing if no data
    if(count==0) return;

    ofstream out(fname);
    for(const auto& sp_name: membr_ptr->species_names){
        fmt::print(out,"{}\t{:.4f}\n",sp_name,around[sp_name]);
    }
    out.close();
}


LipidMembrane::LipidMembrane(System *sys, const std::vector<LipidSpecies> &species, int ngroups)
{
    log = create_logger("membrane");
    system = sys;

    log->info("There are {} lipid species",species.size());
    log->info("Processing lipids...");
    int id = 0;
    for(auto& sp: species){
        vector<Selection> res;
        system->select(sp.whole_str).split_by_residue(res);
        log->info("Lipid {}: {}", sp.name,res.size());

        for(auto& lip: res){
            lipids.emplace_back(lip,sp,id,this);
            ++id;
        }

        species_names.push_back(sp.name);
    }

    // Print statictics
    log->info("Total number of lipids: {}",lipids.size());

    // Create selection for all mid atoms
    all_mid_sel.set_system(*system);
    vector<int> ind;
    ind.reserve(lipids.size());
    for(auto& lip: lipids) ind.push_back(lip.mid_marker_sel.index(0));
    all_mid_sel.modify(ind);

    // Create groups
    groups.reserve(ngroups);
    for(int i=0;i<ngroups;++i) groups.emplace_back(this,i);    
    log->info("{} groups created",ngroups);
}

void LipidMembrane::reset_groups(){
    for(auto& gr: groups){
        gr.reset();
    }
}


void LipidMembrane::compute_properties(float d)
{
    // Set markers for all lipids
    // This unwraps each lipid
    for(auto& l: lipids) l.set_markers();    

    // Get connectivity
    vector<Vector2i> bon;
    vector<float> dist;
    search_contacts(d,all_mid_sel,bon,dist,false,fullPBC);

    // Fill patches with id's and distances
    for(int i=0;i<bon.size();++i){
        int l1 = bon[i](0);
        int l2 = bon[i](1);
        lipids[l1].patch.neib_id.push_back(l2);
        lipids[l1].patch.neib_dist.push_back(dist[l2]);
        lipids[l2].patch.neib_id.push_back(l1);
        lipids[l2].patch.neib_dist.push_back(dist[l1]);
    }

    // Compute local coordinates and approximate normals of patches
    for(size_t i=0; i<lipids.size(); ++i){
        auto& patch = lipids[i].patch;
        // Save central point
        patch.original_center = lipids[i].mid_marker;
        // Set local selection for this lipid
        lipids[i].local_sel = all_mid_sel.select(patch.neib_id);
        lipids[i].local_sel_with_self = lipids[i].local_sel;
        lipids[i].local_sel_with_self.append(all_mid_sel.index(i));

        // Get inertia axes
        Vector3f moments;
        lipids[i].local_sel_with_self.inertia(moments, patch.axes, fullPBC);
        // Transformation matrices
        patch.to_lab = patch.axes.colwise().normalized();
        patch.to_local = patch.to_lab.inverse();
        // Approximate normal with correct orientation
        patch.normal = patch.axes.col(2).normalized();
        float ang = angle_between_vectors(patch.normal, lipids[i].tail_head_vector);
        if(ang > M_PI_2) patch.normal *= -1;
    }

    // Inspect normals and try to fix them if weird orientation is found
    // This is a very important step!
    for(size_t i=0; i<lipids.size(); ++i){
        const Vector3f& normal1 = lipids[i].patch.normal;
        int n_bad = 0;
        Vector3f aver_closest(0,0,0);
        for(int j=0; j<lipids[i].patch.neib_id.size(); ++j){
            int jl = lipids[i].patch.neib_id[j];
            float d = lipids[i].patch.neib_dist[j];
            const Vector3f& normal2 = lipids[jl].patch.normal;
            // We
            if( d<1.0){
                if(angle_between_vectors(normal1,normal2)>M_PI_4 ) ++n_bad;
                aver_closest += normal2;
            }
        }

        if(n_bad>2){
            // Set to average of closest normals
            log->debug("Trying to fix bad normal for lipid {}", i);
            lipids[i].patch.normal = aver_closest.normalized();
        }
    }

    //================
    // Process lipids
    //================
    for(size_t i=0;i<lipids.size();++i){
        LipidMolecule& lip = lipids[i];

        // Sort neib_id for correct index match with selection
        sort(lip.patch.neib_id.begin(),lip.patch.neib_id.end());


        // Create array of local points in local basis
        MatrixXf coord(3,lip.local_sel.size()+1);
        coord.col(0) = Vector3f::Zero(); // Local coord of central point is zero and it goes first
        for(int j=0; j<lip.local_sel.size(); ++j)
            coord.col(j+1) = lip.patch.to_local
                    * system->box(0).shortest_vector(lip.mid_marker,
                                                     lip.local_sel.xyz(j));

        // Fit a quadric surface and find local curvatures
        lip.surf.fit_to_points(coord);
        lip.surf.compute_area();
        lip.surf.compute_curvature();

        // Check direction of the fitted normal
        float norm_ang = angle_between_vectors(lip.patch.to_lab*lip.surf.fitted_normal,
                                               lip.patch.normal);
        if(norm_ang > M_PI_2){
            lip.surf.fitted_normal *= -1.0;
            // Mean curvature have to be inverted as well
            lip.surf.mean_curvature *= -1.0;
            //Gaussian curvature is invariant to the direction of the normal
            // so it should not be flipped!
        }

        lip.mean_curvature = lip.surf.mean_curvature;
        lip.gaussian_curvature = lip.surf.gaussian_curvature;

        // Set neigbours for lipid
        for(int i=0; i<lip.surf.neib_id.size(); ++i){
            // Indexes in lip.surf.neib_id start from 1, because 0 is the central point
            // Thus substract 1 and we get a local selection index
            // which are in turn correspond to lipid indexes in all_mid_sel.
            // The neighbours are NOT arranged in triangulated order!
            int ind = lip.surf.neib_id[i]-1;
            lip.neib.push_back( lip.patch.neib_id[ind] );
        }

        // Set lipid properties
        lip.area = lip.surf.surf_area;
        lip.smoothed_mid_xyz = lip.patch.to_lab*lip.surf.fitted_points.col(0) + lip.mid_marker;
        lip.normal = lip.patch.to_lab*lip.surf.fitted_normal;

        // Tilt
        lip.tilt = rad_to_deg(angle_between_vectors(lip.normal,lip.tail_head_vector));
    } // for lipids

    // Unset markers. This restores all correct atomic coordinates for analysis
    for(auto& lip: lipids) lip.unset_markers();

    // Analysis whith requires restored coordinates of markers
    for(auto& lip: lipids){
        //-------------------------------
        // Dipole
        //-------------------------------
        lip.dipole = lip.whole_sel.dipole(true);
        lip.dipole_proj = lip.dipole.dot(lip.normal); // normal is normalized, no need to devide by norm

        //-----------------------------------
        // Tail properties
        //-----------------------------------
        for(auto& t: lip.tails) t.compute(lip);
    } // Over lipids

    // Process groups
    for(auto& gr: groups){
        gr.process_frame();
    }
}


void LipidMembrane::write_vmd_visualization(const string &path){
    string out1;
    for(const auto& lip: lipids){
        // Print
        const auto& fp = lip.surf.fitted_points;

        Vector3f p1,p2,p3;
        out1 += fmt::format("draw color orange\n");
        // Area vertices in lab coordinates
        out1 += fmt::format("draw color orange\n");
        for(int j=0; j<lip.surf.area_vertexes.size(); ++j){
            int j2 = j+1;
            if(j==lip.surf.area_vertexes.size()-1) j2=0;
            p1 = lip.patch.to_lab *lip.surf.area_vertexes[j] + lip.patch.original_center;
            p2 = lip.patch.to_lab *lip.surf.area_vertexes[j2] + lip.patch.original_center;
            out1 += fmt::format("draw line \"{} {} {}\" \"{} {} {}\" width 2\n",
                       10*p1.x(),10*p1.y(),10*p1.z(),
                       10*p2.x(),10*p2.y(),10*p2.z()
                       );
        }
        // Normal vectors
        // Patch normals
        p1 = lip.patch.to_lab*fp + lip.patch.original_center;
        p2 = p1 + lip.patch.normal*1.0;
        p3 = p1 + lip.patch.normal*1.2;
        out1 += fmt::format("draw color white\n");
        out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2\n",
                   10*p1.x(),10*p1.y(),10*p1.z(),
                   10*p2.x(),10*p2.y(),10*p2.z() );
        out1 += fmt::format("draw cone \"{} {} {}\" \"{} {} {}\" radius 0.3\n",
                   10*p2.x(),10*p2.y(),10*p2.z(),
                   10*p3.x(),10*p3.y(),10*p3.z() );
        // fitted normals
        out1 += fmt::format("draw color iceblue\n");
        p2 = p1 + lip.patch.to_lab*lip.surf.fitted_normal*1.0;
        p3 = p1 + lip.patch.to_lab*lip.surf.fitted_normal*1.3;
        out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2\n",
                   10*p1.x(),10*p1.y(),10*p1.z(),
                   10*p2.x(),10*p2.y(),10*p2.z() );
        out1 += fmt::format("draw cone \"{} {} {}\" \"{} {} {}\" radius 0.3\n",
                   10*p2.x(),10*p2.y(),10*p2.z(),
                   10*p3.x(),10*p3.y(),10*p3.z() );

        // Smoothed atom
        out1 += fmt::format("draw sphere \"{} {} {}\" radius 2\n",
                   10*p1.x(),10*p1.y(),10*p1.z() );
    }

    // Output area and normals plots
    ofstream out_all(fmt::format("{}/areas_all.tcl",path));
    out_all << out1;
    out_all.close();

    for(int i=0; i<lipids.size(); ++i){
        all_mid_sel.beta(i) = 10*lipids[i].mean_curvature;
    }
    all_mid_sel.write(fmt::format("{}/areas_all.pdb",path));
}


void LipidMembrane::compute_averages()
{
    // Run post-processing for all groups
    for(auto& gr: groups){
        gr.post_process();
    }
}

void LipidMembrane::write_averages(string path)
{
    string s;
    s += "Run summary:\n";
    s += fmt::format("Lipid species ({}): {}\n",species_names.size(),fmt::join(species_names," "));
    for(auto& gr: groups){
        s += gr.summary();
    }

    // Print summary
    cout << s << endl;

    // Save summary to file
    ofstream out(path+"/summary.dat");
    out << s;
    out.close();

    // Write files for species properties
    for(int g=0;g<groups.size();++g){
        groups[g].save_properties_table_to_file(fmt::format("{}/gr{}_properties.dat",path,g));

        for(auto& sp: groups[g].species_properties){
            if(sp.second.count>0){
                string file_prefix(fmt::format("{}/gr{}_{}_",path,g,sp.first));
                // Area
                sp.second.area_hist.save_to_file(file_prefix+"area.dat");
                // Tilt
                sp.second.tilt_hist.save_to_file(file_prefix+"tilt.dat");
                // Curvature
                sp.second.mean_curv_hist.save_to_file(file_prefix+"mean_curv.dat");
                sp.second.gauss_curv_hist.save_to_file(file_prefix+"gauss_curv.dat");
                // Order
                sp.second.save_order_to_file(file_prefix+"order.dat");
                // Around
                sp.second.save_around_to_file(file_prefix+"around.dat");
            }
        }
    }
}

LipidGroup::LipidGroup(LipidMembrane *ptr, int id){
    gr_id = id;
    membr_ptr = ptr;
    num_lipids = 0;
    num_frames = 0;
    trans_dihedrals_ratio.fill(0.0);
    // Initialize species_properties
    for(const auto& sp_name: membr_ptr->species_names){
        species_properties.emplace(sp_name,PerSpeciesProperties(membr_ptr));
    }
}

void LipidGroup::process_frame()
{
    // Cycle over lipids in this group
    for(int id: lip_ids){
        auto name = membr_ptr->lipids[id].name;
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
        s += fmt::format("\tTr.Dih.:\t{} +/- {}\n", trans_dihedrals_ratio[0],trans_dihedrals_ratio[1]);
        s += "\tLipid species:\n";
        for(auto& sp: membr_ptr->species_names){
            s += fmt::format("\t{}:\n", sp);
            s += species_properties.at(sp).summary();
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
    for(auto& sp: membr_ptr->species_names){
        s += fmt::format("{}", sp);
        const auto& prop = species_properties.at(sp);
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


void QuadSurface::fit_to_points(const MatrixXf &coord){
    int N = coord.cols();
    // We fit with polynomial fit = A*x^2 + B*y^2 + C*xy + D*x + E*y + F
    // Thus we need a linear system of size 6
    Matrix<float,6,6> m;
    Matrix<float,6,1> rhs; // Right hand side and result
    m.fill(0.0);
    rhs.fill(0.0);

    Vector<float,6> powers;
    powers(5) = 1.0; //free term, the same everywhere
    for(int j=0;j<N;++j){ // over points
        powers(0) = coord(0,j)*coord(0,j); //xx
        powers(1) = coord(1,j)*coord(1,j); //yy
        powers(2) = coord(0,j)*coord(1,j); //xy
        powers(3) = coord(0,j); //x
        powers(4) = coord(1,j); //y
        // Fill the matrix
        m.noalias() += powers * powers.transpose();
        // rhs
        rhs += powers*coord(2,j);
    }

    // Now solve
    quad_coefs = m.ldlt().solve(rhs);

    // Compute RMS of fitting and fitted points
    fitted_points = coord; // Set to parent points initially
    fit_rms = 0.0;
    for(int j=0;j<N;++j){
        float fitZ = evalZ(coord(0,j),coord(1,j));
        fit_rms += pow(coord(2,j)-fitZ,2);
        // Adjust Z of fitted point
        fitted_points(2,j) = fitZ;
    }
    fit_rms = sqrt(fit_rms/float(N));
}

void QuadSurface::compute_area(){
    // The center of the cell is assumed to be at {0,0,0}
    // First point is assumed to be the central one and thus not used
    using namespace voro;

    voronoicell_neighbor cell;
    // Initialize with large volume by default in XY and size 1 in Z
    cell.init(-10,10,-10,10,-0.5,0.5);
    for(int i=1; i<fitted_points.cols(); ++i){ // From 1, central point is 0 and ignored
        cell.nplane(fitted_points(0,i),fitted_points(1,i),0.0,i); // Pass i as neigbour pid
    }

    vector<int> neib_list;
    vector<int> face_vert;
    vector<double> vert_coords;

    //cell.draw_gnuplot(0,0,0,"cell.dat");

    cell.neighbors(neib_list);

    cell.face_vertices(face_vert);
    // In face_vert the first entry is a number k corresponding to the number
    // of vertices making up a face,
    // this is followed by k additional entries describing which vertices make up this face.
    // For example, (3, 16, 20, 13) is a triangular face linking vertices 16, 20, and 13

    cell.vertices(0,0,0,vert_coords);
    //vert_coords this corresponds to a vector of triplets (x, y, z)
    // describing the position of each vertex.

    // Collect verteces which make up cell area in XY plane
    int j = 0;
    int vert_ind;
    float x,y;
    area_vertexes.clear();
    area_vertexes.reserve(cell.number_of_faces()-2); // 2 faces are top and bottom walls

    // We need to take only face aganst the wall, which containes all vertices
    // which we need to area calculation
    for(int i=0;i<neib_list.size();++i){
        if(neib_list[i]<0){ // Only take the wall faces
            // Cycle over vertices of this face
            for(int ind=0; ind<face_vert[j]; ++ind){
                vert_ind = 3*face_vert[j+1+ind];
                x = vert_coords[vert_ind];
                y = vert_coords[vert_ind+1];
                // We are not interested in Z component
                area_vertexes.push_back(Vector3f(x,y,0.0));
            }
            break; // We are done
        }
        j += face_vert[j]+1; // Next entry in face_vert
    }

    // In plane area. Since Z size of the cell is 1 volume=area
    in_plane_area = cell.volume();

    // Compute area on surface
    // Project all vertexes to the surface
    surf_area = 0.0;
    for(Vector3f& v: area_vertexes) project_point_to_surface(v);
    // Sum up areas of triangles. Points are not dupicated.
    for(int i=0; i<area_vertexes.size(); ++i){
        int ii = (i<area_vertexes.size()-1) ? i+1 : 0;
        surf_area += 0.5*(area_vertexes[i]-fitted_points.col(0))
                .cross(area_vertexes[ii]-fitted_points.col(0))
                .norm();
    }

    // Set list of neigbours
    // Filter out negatives
    for(int id: neib_list){
        if(id>=0) neib_id.push_back(id);
    }
}

void QuadSurface::compute_curvature(){
    /* Compute the curvatures

            First fundamental form:  I = E du^2 + 2F du dv + G dv^2
            E= r_u dot r_u, F= r_u dot r_v, G= r_v dot r_v

            For us parametric variables (u,v) are just (x,y) in local space.
            Derivatives:
                r_u = {1, 0, 2Ax+Cy+D}
                r_v = {0, 1, 2By+Cx+E}

            In central point x=0, y=0 so:
                r_u={1,0,D}
                r_v={0,1,E}

            Thus: E_ =1+D^2; F_ = D*E; G_ = 1+E^2;

            Second fundamental form: II = L du2 + 2M du dv + N dv2
            L = r_uu dot n, M = r_uv dot n, N = r_vv dot n

            Normal is just  n = {0, 0, 1}
            Derivatives:
                r_uu = {0, 0, 2A}
                r_uv = {0 ,0, C}
                r_vv = {0, 0, 2B}

            Thus: L_ = 2A; M_ = C; N_ = 2B;
            */

    float E_ = 1.0+D()*D();
    float F_ = D()*E();
    float G_ = 1.0+E()*E();

    float L_ = 2.0*A();
    float M_ = C();
    float N_ = 2.0*B();

    //Curvatures:
    gaussian_curvature = (L_*N_-M_*M_)/(E_*G_-F_*F_);
    mean_curvature = 0.5*(E_*N_-2.0*F_*M_+G_*L_)/(E_*G_-F_*F_);

    // Compute normal of the fitted surface at central point
    // dx = 2Ax+Cy+D
    // dy = 2By+Cx+E
    // dz = -1
    // Since we are evaluating at x=y=0:
    // norm = {D,E,1}
    fitted_normal = Vector3f(D(),E(),-1.0).normalized();
    // Orientation of the normal could be wrong!
    // Have to be flipped according to lipid orientation later
}

/*
// Triangulated surface output
string out2;
Vector3f p1,p2,p3,n1,n2,n3;
for(int i=0; i<lipids.size(); ++i){
    for(int j=0; j<lipids[i].neib.size(); ++j){
        int ind2 = lipids[i].neib[j];
        int k = (j<lipids[i].neib.size()-1) ? j+1 : 0;
        int ind3 = lipids[i].neib[k];
        p1 = all_mid_sel.xyz(i);
        p2 = all_mid_sel.box().closest_image(all_mid_sel.xyz(ind2),p1);
        p3 = all_mid_sel.box().closest_image(all_mid_sel.xyz(ind3),p1);
        n1 = lipids[i].normal;
        n2 = lipids[j].normal;
        n3 = lipids[k].normal;

        out2 += fmt::format("draw sphere \"{} {} {}\" radius 2\n",
                   10*p1.x(),10*p1.y(),10*p1.z() );

        out2 += fmt::format("draw tricolor "
                            "\"{} {} {}\" \"{} {} {}\" \"{} {} {}\" "
                            "\"{} {} {}\" \"{} {} {}\" \"{} {} {}\" "
                            "{} {} {} \n",
                            10*p1.x(),10*p1.y(),10*p1.z(),
                            10*p2.x(),10*p2.y(),10*p2.z(),
                            10*p3.x(),10*p3.y(),10*p3.z(),

                            n1.x(),n1.y(),n1.z(),
                            n2.x(),n2.y(),n2.z(),
                            n3.x(),n3.y(),n3.z(),

                            i, j, k
                            );
    }
}
ofstream outf2(fmt::format("vis/surf.tcl"));
outf2 << out2;
outf2.close();
*/
