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

#include "pteros/extras/membrane/lipid_membrane.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>
#include <unordered_set>
#include "voro++.hh"
#include <omp.h>
#include <filesystem>
#include "fmt/os.h"
#include <ranges>
#include <set>
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

struct InterpData {
    int lip_ind; // Lipid index
    float surf_dist; // Distance from lipid surf marker
    float norm_dist; // Distance from lipid normal
    float depth; // Depth along the normal
};


float gaussian(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;
    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}


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

LipidMembrane::LipidMembrane(const Selection &input_sel,
                             int ngroups,
                             const vector<LipidSpecies>& sp_list,
                             const Selection &incl):
    input_sel_ptr(&input_sel),
    inclusion(incl)
{
    log = spdlog::get("membrane");
    if(!log) log = create_logger("membrane");

    // Add all lipid species provided
    species.reserve(sp_list.size()); // Reserve space to avoid reallocation and ensure correct pointers!
    for(auto& sp: sp_list) register_lipid_species(sp);

    // Create groups (they see all registered species)
    groups.reserve(ngroups);
    for(int i=0;i<ngroups;++i) groups.emplace_back(this,i);
    log->info("{} groups created",ngroups);

    // Print statictics
    log->info("Total number of lipids: {}",lipids.size());

    // Create internal system for surf markers
    for(size_t i=0; i<lipids.size(); ++i){
        Atom at;
        at.name = "SRF";
        at.resname = "SRF";
        at.resid = lipids[i].whole_sel.resid(0);
        at.mass = 1.0;
        surf_sys.append(at,Vector3f::Zero());
    }

    // Create selection for all surf marker
    all_surf_sel = surf_sys.select_all();

    // If asked for per-carbon normals make a selection of all tail carbons
    if(options.per_carbon_normals){
        all_tails_sel.set_system(*(input_sel_ptr->get_system()));
        for(auto& lip: lipids){
            for(auto& t: lip.tails){
                t.first_global_index = all_tails_sel.size();
                // Add indexes of this tail to global vector
                for(int offset: t.get_descr().c_offsets){
                    // We use unchecked variant to create selection index without sorting
                    // exactly in the order of addition
                    all_tails_sel.append_unchecked(lip.whole_sel.index(offset));
                }
            }
        }

        log->info("Per-atom normals requested for {} tail carbons",all_tails_sel.size());
    }
}


std::vector<Selection> LipidMembrane::get_domains(const System &sys, const std::vector<LipidSpecies> &sp_list, float d)
{
    std::string domain_names_str;
    Selection domains_sel;
    std::vector<Selection> domains;

    auto dom_log = create_logger("get_domains");

    for(auto const& sp: sp_list){
        // Extract tail carbon names and add them to domain_names_str
        string tnames = "";
        for(auto const& tstr: sp.tails_descr_strings){
            size_t tok_start = 0;
            size_t cur = 0;
            while(cur<tstr.size()){
                if(tstr[cur]=='-' || tstr[cur]=='='){
                    tnames += tstr.substr(tok_start,cur-tok_start)+" ";
                    tok_start=cur+1;
                }
                ++cur;
            }
            tnames += tstr.substr(tok_start);
        }
        if(domain_names_str!="") domain_names_str += " or ";
        domain_names_str += fmt::format("({} and name {})",sp.whole_sel_str,tnames);
    }

    // Find domains
    // Create domains selection
    //log->info("{}",domain_names_str);

    domains_sel = sys.select(domain_names_str);
    dom_log->info("Will find domains using {} tail atoms",domains_sel.size());
    domains_sel.split_by_connectivity(d,domains,true);
    dom_log->info("There are {} domains",domains.size());

    vector<Selection> out;
    out.reserve(domains.size());
    for(auto const& d: domains){
        string s="";
        for(auto const& r: d.get_resindex(true)) s+=fmt::format("{} ",r);
        out.emplace_back(sys,fmt::format("by residue (resindex {})",s));
    }

    return out;
}


void LipidMembrane::register_lipid_species(const LipidSpecies &sp)
{
    log->info("Adding lipid species {}",sp.name);
    species.push_back(sp);

    vector<Selection> res;
    input_sel_ptr->select(sp.whole_sel_str).split_by_residue(res);
    if(res.size()==0){
        log->warn("Skipping {} - there are no such lipids in the system!",sp.name);
        return;
    }
    log->info("There are {} {} lipids in the system", res.size(),sp.name);

    // Add lipids
    int id = lipids.size();
    for(size_t i=0; i<res.size(); ++i){
        // For the first lipid of this kind initialize
        // structure-dependent data in LipidSpecies
        if(i==0) species.back().init(res[i]);
        // Add lipid
        lipids.emplace_back(res[i],&species.back(),id,this);
        ++id;
    }
}


void LipidMembrane::reset_groups(){
    for(auto& gr: groups){
        gr.reset();
    }
}


void LipidMembrane::get_initial_normals(){
    // First pass: average of head-tail vector over a patch (including self)
    for(size_t i=0; i<lipids.size(); ++i){
        auto& lip = lipids[i];
        lip.patch.normal = lip.tail_head_vector;
        for(size_t i=0; i<lip.patch.neib_id.size(); ++i){
            size_t ind = lip.patch.neib_id[i];
            //float w = lip.inclusion_neib.empty() ? gaussian(patch.neib_dist[i],0,2.0) : 1.0;
            lip.patch.normal += lipids[ind].tail_head_vector; // * w;

        }
        lip.patch.normal.normalize();
    }

    // Second pass: average of normals from the first pass over a patch (including self)
    for(size_t i=0; i<lipids.size(); ++i){
        auto& lip = lipids[i];
        for(size_t i=0; i<lip.patch.neib_id.size(); ++i){
            size_t ind = lip.patch.neib_id[i];
            //float w = lip.inclusion_neib.empty() ? gaussian(patch.neib_dist[i],0,2.0) : 1.0;
            lip.patch.normal += lipids[ind].patch.normal; // * w;

        }
        lip.patch.normal.normalize();

        // Correct orientation
        float ang = angle_between_vectors(lip.patch.normal, lip.tail_head_vector);
        if(ang > M_PI_2) lip.patch.normal *= -1;
    }
}

void LipidMembrane::create_lipid_patches(float d){
    vector<Vector2i> bon;
    vector<float> dist;
    search_contacts(d,all_surf_sel,bon,dist,false,fullPBC);

    for(auto& lip: lipids) {
        // Clear patches for all lipids
        lip.patch.neib_id.clear();
    }

    // Fill patches with id's
    for(auto const& b: bon){
        int l1 = b(0);
        int l2 = b(1);
        lipids[l1].patch.neib_id.push_back(l2);
        lipids[l2].patch.neib_id.push_back(l1);
    }

    // Sort neib_id of patches and remove duplicates for correct index match with all_surf_sel selection
    for(auto& lip: lipids) {
        lip.patch.sort_and_remove_duplicates();
    }
}


void LipidMembrane::add_inclusions(float incl_d){
    //fmt::print("Adding inclusion {}\n",inclusion.size());
    inclusion.apply(); // In case if it is coord-dependent
    vector<Vector2i> bon;
    vector<float> dist;
    search_contacts(incl_d,all_surf_sel,inclusion,bon,dist,false,fullPBC);

    // For each lipid add contacting inclusion atoms
    for(size_t i=0;i<bon.size();++i){
        int l = bon[i](0);
        int in = bon[i](1);
        lipids[l].inclusion_neib.push_back(in);
    }
}

void LipidMembrane::smoothed_markers_to_surf_coord(int ind){
    auto const& box = input_sel_ptr->box();
    auto& lip = lipids[ind];
    size_t n = lip.patch.neib_id.size();
    // Allocate array of local points in local basis
    lip.surf.lip_coord.resize(3,n+1);
    // Central point is zero and it goes first
    lip.surf.lip_coord.col(0) = Vector3f::Zero();
    for(size_t j=0; j<n; ++j){
        int id = lip.patch.neib_id[j];
        lip.surf.lip_coord.col(j+1) = box.shortest_vector(lip.smoothed_surf_marker,
                                                          lipids[id].smoothed_surf_marker);
    }
}


void LipidMembrane::inclusion_coord_to_surf_coord(int ind){
    auto const& box = input_sel_ptr->box();
    auto& lip = lipids[ind];
    lip.surf.inclusion_coord.resize(3,lip.inclusion_neib.size());
    for(size_t i=0; i<lip.inclusion_neib.size(); ++i){
        lip.surf.inclusion_coord.col(i) = box.shortest_vector(
                lip.smoothed_surf_marker,
                inclusion.xyz(lip.inclusion_neib[i])
              );

    }
}

void LipidMembrane::do_smoothing(int i, Vector3f_const_ref normal){
    // Create array of local points
    smoothed_markers_to_surf_coord(i);
    // Compute fitting polynomial and update fitted normal and central point
    lipids[i].surf.fit_quad(normal);
    // Set smoothed point in lab coordinates
    lipids[i].smoothed_surf_marker += lipids[i].surf.fitted_central_point;
}

void LipidMembrane::lip_neib_from_voronoi(int l){
    auto& lip = lipids[l];
    lip.neib.clear();
    for(size_t i=0; i<lip.surf.neib_id.size(); ++i){
        // Indexes in lip.surf.neib_id start from 1, because 0 is the central point
        // Thus substract 1 and we get a local selection index
        // which are in turn correspond to lipid indexes in all_surf_sel.
        // The neighbours are NOT arranged in triangulated order!
        int ind = lip.surf.neib_id[i]-1;
        lip.neib.push_back( lip.patch.neib_id[ind] );
    }
}

//---------------------------------------------------------------------------------

void LipidMembrane::compute_properties(float d, float incl_d, OrderType order_type)
{
    auto const& box = input_sel_ptr->box();

    // Update box in surf_sys
    surf_sys.box() = box;

    // Set markers for all lipids
    // This unwraps each lipid    
    for(size_t i=0; i<lipids.size(); ++i){
        lipids[i].set_markers();
        all_surf_sel.xyz(i) = lipids[i].surf_marker;
    }

    // Populate lipid patches
    create_lipid_patches(d);

    // Guess initial normals
    get_initial_normals();

    //#pragma omp parallel for if (lipids.size() >= 100)

    //=====================
    // Iterative smoothing
    //=====================

    // Iterations update the smoothed_surf_marker,
    // so set it to initial surf market for the first iteration
    for(auto& lip: lipids){
        lip.smoothed_surf_marker = lip.surf_marker;
    }

    // Do initial smoothing with patch normals
    for(size_t i=0;i<lipids.size();++i){
        do_smoothing(i,lipids[i].patch.normal);
    }

    // Iterations    
    size_t iter = 0;
    float rmsd, prev_rmsd=1e6;    
    vector<QuadSurface> prev_state(lipids.size());

    do{
        rmsd = 0.0;
        for(size_t i=0;i<lipids.size();++i){
            auto& lip = lipids[i];
            prev_state[i] = lip.surf;
            Vector3f prev_marker = lip.smoothed_surf_marker;
            // Use existing fitted normal
            do_smoothing(i,lip.surf.fitted_normal);
            rmsd += (lip.smoothed_surf_marker-prev_marker).squaredNorm();
        }
        ++iter;
        rmsd = sqrt(rmsd/lipids.size());

        log->debug("Iter={} rmsd={}",iter,rmsd);

        if(rmsd>prev_rmsd){
            // Oversmoothing starts, revert to saved markers and exit the loop
            for(size_t i=0;i<lipids.size();++i){
                lipids[i].surf = prev_state[i];
                lipids[i].smoothed_surf_marker = lipids[i].surf.fitted_central_point;
            }
            break;
        } else {
            prev_rmsd = rmsd;
        }

    } while(rmsd>options.smoothing_tol && iter<options.smoothing_maxiter);

    //=====================
    // Lipid properties
    //=====================

    // If inclusions are present compute contacts with them
    if(inclusion.size()>0){
        add_inclusions(incl_d);
    }

    for(size_t i=0;i<lipids.size();++i){
        auto& lip = lipids[i];

        // If inclusion is present bring neibouring inclusion atoms
        // to local basis as well
        if(!lip.inclusion_neib.empty()){
            inclusion_coord_to_surf_coord(i);
        }

        lip.surf.compute_voronoi(options.inclusion_h_cutoff);

        // area
        lip.area = lip.surf.surf_area;

        // If area is zero then the lipid is invalid, so stop here
        if(lip.area==0) continue;

        //lip.surf.compute_curvature_and_normal();

        // Check direction of the fitted normal
        float norm_ang = angle_between_vectors(lip.surf.fitted_normal,lip.patch.normal);
        if(norm_ang > M_PI_2){
            lip.surf.fitted_normal *= -1.0;
            // Mean curvature have to be inverted as well
            lip.surf.mean_curvature *= -1.0;
            //Gaussian curvature is invariant to the direction of the normal
            // so it should not be flipped!
        }
        // Curvatures
        lip.mean_curvature = lip.surf.mean_curvature;
        lip.gaussian_curvature = lip.surf.gaussian_curvature;
        // Normal
        lip.normal = lip.surf.fitted_normal;
        // Tilt
        float a = angle_between_vectors(lip.normal,lip.tail_head_vector);
        lip.tilt = rad_to_deg(a);
        // Monolayer thickness
        lip.mono_thickness = lip.tail_head_vector.norm()*cos(a);

        // Set neigbours for lipid
        lip_neib_from_voronoi(i);

        // Coordination number
        lip.coord_number = lip.neib.size();
    } // for lipids

    // If asked make normal interpolation for all tail atoms used for order calculation
    MatrixXf c_normals;
    if(options.per_carbon_normals){
        vector<InterpolatedPoint> interp;        
        get_interpolation(all_tails_sel,interp,c_normals);
    }

    for(auto& lip: lipids){
        // Tail properties
        for(auto& t: lip.tails){
            if(options.per_carbon_normals){
                // We need to check the normal orientation because the ends of the
                // tails may have flipped normals due to the influence of the opposite monolayer
                // For this we need to pass the lipid normal as a reference
                t.compute_order_and_dihedrals(lip.whole_sel,
                                              c_normals,
                                              lip.normal,
                                              order_type);
            } else {
                t.compute_order_and_dihedrals(lip.whole_sel,
                                              lip.normal,
                                              lip.normal,
                                              order_type);
            }
        }
    }

    // Process groups
    for(auto& gr: groups){
        gr.process_frame();
    }
}


void LipidMembrane::get_interpolation(const Selection &points,
                                      vector<InterpolatedPoint>& res,
                                      MatrixXf& normals, float d)
{
    // Resize coeffs correctly
    res.resize(points.size());
    normals.resize(3,points.size());

    // Find distances between provided points and lipid markers
    // get local indexes
    vector<Vector2i> bon;
    vector<float> dist;
    search_contacts(d, points, all_surf_sel,
                    bon,dist,false,fullPBC);

    // Convert to i->(j,k,l..) keeping distances as well
    // Indexed by point, contains lipid indexes
    //                 ind  dist
    vector<vector<InterpData>> con(points.size());
    for(size_t i=0;i<bon.size();++i){
        //                        lipid      dist     not used yet
        con[bon[i](0)].push_back({bon[i](1), dist[i], 0.0, 0.0});
    }

    for(int p=0; p<points.size(); ++p){

        if(con[p].size()<1){
            // This point is too far from the markers
            // Perform an exhaustive search and take only one nearest lipid
            float mind = std::numeric_limits<float>::max();
            int mini = -1;
            auto const& box = all_surf_sel.box();
            for(int i=0; i<all_surf_sel.size(); ++i){
                float d = box.distance(points.xyz(p),all_surf_sel.xyz(i));
                if(d<mind){
                    mini = i;
                    mind = d;
                }
            }
            con[p].push_back({mini, mind, 0.0, 0.0});
        }

        // Do normal filtering if more than one neighbour
        if(con[p].size()>1){
            // Find the closest lipid
            int closest_ind = -1;
            float min_dist = std::numeric_limits<float>::max();
            for(size_t i=0;i<con[p].size();++i){
                if(con[p][i].surf_dist<min_dist){
                    min_dist = con[p][i].surf_dist;
                    closest_ind = i;
                }
            }

            // Keep only those lipids which have the same normal orientation as the closest one
            // (filter the accidental inclusion of the opposite monolayer)
            auto const& n0 = lipids[con[p][closest_ind].lip_ind].normal;
            for(auto& el: con[p]){
                if(angle_between_vectors(lipids[el.lip_ind].normal,n0)>M_PI_2){
                    el.lip_ind = -1; // Indicates that this lipid is not valid
                }
            }
        }

        // For all valid lipids compute distance from normal and depth
        auto const& box = input_sel_ptr->box();
        for(auto& el: con[p]){
            if(el.lip_ind>=0){
                auto const& x1 = all_surf_sel.xyz(el.lip_ind);
                auto const& x2 = x1+lipids[el.lip_ind].normal;
                auto const& x0 = points.xyz(p);
                // Formula for distance from: https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
                // no need to divide by /(x2-x1).norm() because lipid normal is already normalized
                el.norm_dist = (box.shortest_vector(x1,x0).cross(box.shortest_vector(x2,x0))).norm();
                // Formula for scalar projection from: https://en.wikipedia.org/wiki/Dot_product#Scalar_projection_and_first_properties
                el.depth = box.shortest_vector(x1,x0).dot(lipids[el.lip_ind].normal);
            }
        }

        float sumw = 0.0;
        float sumw_c = 0.0;
        float mean_c = 0.0;
        for(auto& el: con[p]){
            if(el.lip_ind>=0){
                // Small gaussian width to include only very close normals
                float w = gaussian(el.norm_dist,0.0,0.5);
                // For curvature use broader gaussian for better averaging
                float w1 = gaussian(el.surf_dist,0.0,2.5);
                sumw += w;
                sumw_c += w1;
                res[p].neib_lipids.push_back(el.lip_ind);
                res[p].weights.push_back(w);
                res[p].normal += lipids[el.lip_ind].normal * w;
                res[p].mean_depth += w * el.depth;
                mean_c += w1 * lipids[el.lip_ind].mean_curvature;
            }
        }

        // Normalize
        for(auto& v: res[p].weights) v /= sumw;
        res[p].normal /= sumw;
        res[p].mean_depth /= sumw;
        mean_c /= sumw_c;

        normals.col(p) = res[p].normal;

        if(mean_c>0)
            res[p].mean_curvature = mean_c/ (1.0-mean_c * res[p].mean_depth);
        else
            res[p].mean_curvature = mean_c/ (1.0+mean_c * res[p].mean_depth);
    } // over points
}

MatrixXf LipidMembrane::get_average_curvatures(int lipid, int n_shells)
{
    MatrixXf m(n_shells,2);
    m.fill(0.0);

    vector<int> neib_n;
    neib_n.push_back(lipid); //self

    unordered_set<int> s;
    // Compute up to order n_shells neighbours
    for(int n=0; n<n_shells; ++n){
        // Compute averages
        for(int el: neib_n){
            m(n,0) += lipids[el].mean_curvature;
            m(n,1) += lipids[el].gaussian_curvature;
        }
        m.row(n) /= float(neib_n.size());

        // Update neib2
        for(int n1: neib_n){
            for(int n2: lipids[n1].neib) s.insert(n2);
        }
        std::copy(s.begin(),s.end(),back_inserter(neib_n));
    }

    return m;
}


void LipidMembrane::compute_triangulation()
{
    const auto& box = lipids[0].whole_sel.box();

    // Convert neibour lists to sets for fast search
    vector<unordered_set<int>> neib(lipids.size());
    for(size_t i=0; i<lipids.size(); ++i){
        for(auto el: lipids[i].neib) neib[i].insert(el);
    }

    vector<Vector3i> triangles;

    // Now perform search
    //for(int i1=0; i1<lipids.size(); ++i1){
    for(int i1=40; i1<41; ++i1){
        for(int i2: neib[i1]){ // Loop over neibours of i1 and get i2
            // We only work with i1<i2 to avoid double counting
            //if(i1>i2) continue;
            // Loop over neibours of i2 and get i3
            for(int i3: neib[i2]){
                // We only work with i2<i3 to avoid double counting
                //if(i2>i3) continue;
                // Look if i3 is within neigbours of i1
                auto it = neib[i1].find(i3);
                if(it!=neib[i1].end()){
                    // Found, which means that i1-i2-i3 is a triangle
                    // See normal orientation
                    Vector3f v1 = box.shortest_vector(lipids[i1].smoothed_surf_marker,lipids[i2].smoothed_surf_marker);
                    Vector3f v2 = box.shortest_vector(lipids[i1].smoothed_surf_marker,lipids[i3].smoothed_surf_marker);
                    Vector3f n = v1.cross(v2);
                    if(angle_between_vectors(n,lipids[i1].normal)<M_PI_2){
                        triangles.push_back({i1,i2,i3});
                    } else {
                        triangles.push_back({i3,i2,i1});
                    }
                }
            }
        }

    }

    /*
    // Output an obj file
    string s;
    // Vertexes
    for(auto l: lipids){
        s+=fmt::format("v {} {} {}\n",l.mid_marker.x(),l.mid_marker.y(),l.mid_marker.z());
    }
    // Normals
    for(auto l: lipids){
        s+=fmt::format("vn {} {} {}\n",l.normal.x(),l.normal.y(),l.normal.z());
    }
    // Triangular faces
    const auto& b = lipids[0].whole_sel.box();
    for(auto t: triangles){
        // Only output face if vertexes are not pbc-wrapped
        Vector3f d1 = lipids[t(0)].mid_marker-lipids[t(1)].mid_marker;
        Vector3f d2 = lipids[t(0)].mid_marker-lipids[t(2)].mid_marker;
        if(   d1(0)<0.5*b.extent(0) && d1(1)<0.5*b.extent(1) && d1(2)<0.5*b.extent(2)
           && d2(0)<0.5*b.extent(0) && d2(1)<0.5*b.extent(1) && d2(2)<0.5*b.extent(2)
          ){
            s+=fmt::format("f {}//{} {}//{} {}//{}\n",t(0)+1,t(0)+1,t(1)+1,t(1)+1,t(2)+1,t(2)+1);
        }
    }
    */

    // Output VMD surface
    // Prepare VMD indexed colors
    // We use linear interpolation between start, mid and end.
    Vector3f c1{0,0,1}; //Blue
    Vector3f c2{1,1,1}; //White
    Vector3f c3{1,0,0}; //Red
    int n_colors = 104;
    MatrixXf colors(3,n_colors);
    for(int i=0; i<n_colors/2; ++i){
        colors.col(i) = c1 + i*(c2-c1)/float(n_colors/2-1);
    }
    for(int i=0; i<n_colors/2; ++i){
        colors.col(i+n_colors/2) = c2 + i*(c3-c2)/float(n_colors/2-1);
    }

    vector<MatrixXf> curv(lipids.size());
    for(size_t i=0; i<lipids.size(); ++i){
        curv[i] = get_average_curvatures(i,5);
    }

    for(int smooth_level=0;smooth_level<1;++smooth_level){
        // Now assign color to each lipid according to mean curvature
        float min_c=1e6;
        float max_c=-1e6;
        for(size_t i=0; i<lipids.size(); ++i){
            float c = curv[i](smooth_level,0);
            if(c<min_c) min_c = c;
            if(c>max_c) max_c = c;
        }
        fmt::print("min_c={} max_c={}\n",min_c,max_c);

        // curv           x     x=103*(curv-min_c)/(max_c-min_c)
        // max_c-min_c    104
        VectorXi color_ind(lipids.size());
        for(size_t i=0; i<lipids.size(); ++i){
            color_ind[i] = 1057 - n_colors +
                    round((n_colors-1)*(curv[i](smooth_level,0)-min_c)/(max_c-min_c));
        }

        // Update VMD color definitions
        string s;
        for(int i=0; i<n_colors; ++i){
            s+= fmt::format("color change rgb {} {} {} {}\n",
                            1057-n_colors+i,
                            colors(0,i),
                            colors(1,i),
                            colors(2,i));
        }


        // Output all neib edges

        s+="draw materials on\n";
        s+="draw material Diffuse\n";

        for(size_t i=0; i<lipids.size(); ++i){
            Vector3f p1 = lipids[i].smoothed_surf_marker;
            s+=fmt::format("draw color {}\n",color_ind[i]);
            s+=fmt::format("draw sphere \"{} {} {}\" radius 1.3 resolution 12\n",
                           p1(0)*10,p1(1)*10,p1(2)*10);

            s+=fmt::format("draw color violet\n");
            for(int nb: lipids[i].neib){
                Vector3f p2 = box.closest_image(lipids[nb].smoothed_surf_marker,p1);
                s+=fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2\n",
                               p1(0)*10,p1(1)*10,p1(2)*10,
                               p2(0)*10,p2(1)*10,p2(2)*10
                              );
            }

        }

        // Output colored triangles
        for(auto t: triangles){
            Vector3f p1 = lipids[t(0)].smoothed_surf_marker;
            Vector3f p2 = lipids[t(1)].smoothed_surf_marker;
            Vector3f p3 = lipids[t(2)].smoothed_surf_marker;
            p2 = box.closest_image(p2,p1);
            p3 = box.closest_image(p3,p1);
            p1*=10.0;
            p2*=10.0;
            p3*=10.0;
            Vector3f n1 = lipids[t(0)].normal.normalized();
            Vector3f n2 = lipids[t(1)].normal.normalized();
            Vector3f n3 = lipids[t(2)].normal.normalized();
            int c1 = color_ind[t(0)];
            int c2 = color_ind[t(1)];
            int c3 = color_ind[t(2)];

            s += fmt::format("draw tricolor "
                    "\"{} {} {}\" " //p1
                    "\"{} {} {}\" " //p2
                    "\"{} {} {}\" " //p3
                    "\"{} {} {}\" " //n1
                    "\"{} {} {}\" " //n2
                    "\"{} {} {}\" " //n3
                    "{} {} {}\n", //c1 c2 c3
                    p1(0),p1(1),p1(2),
                    p2(0),p2(1),p2(2),
                    p3(0),p3(1),p3(2),
                    n1(0),n1(1),n1(2),
                    n2(0),n2(1),n2(2),
                    n3(0),n3(1),n3(2),
                    c1,c2,c3);

            /*
            s += fmt::format("draw triangle "
                    "\"{} {} {}\" " //p1
                    "\"{} {} {}\" " //p2
                    "\"{} {} {}\"\n", //p3
                     p1(0),p1(1),p1(2),
                     p2(0),p2(1),p2(2),
                     p3(0),p3(1),p3(2)
                    );
                    */
        }

        auto f = fmt::output_file(fmt::format("triangulated_smooth_level_{}.tcl",smooth_level));
        f.print("{}",s);
        f.close();
    } //smooth level
    //exit(1);
}


void LipidMembrane::write_vmd_visualization(const string &path){
    string out1;

    for(const auto& lip: lipids){
        // Skip invalid lipids
        if(lip.area==0) continue;

        Vector3f p1,p2,p3;

        // Voronoi vertices in lab coordinates

        out1+="draw materials on\n";
        out1+="draw material AOEdgy\n";
        out1 += fmt::format("draw color orange\n");
        for(size_t j=0; j<lip.surf.area_vertexes.size(); ++j){
            int j2 = j+1;
            if(j==lip.surf.area_vertexes.size()-1) j2=0;
            p1 = lip.surf.area_vertexes[j]  + lip.smoothed_surf_marker;
            p2 = lip.surf.area_vertexes[j2] + lip.smoothed_surf_marker;
            // Draw small sphere at cylinders intersection
            out1 += fmt::format("draw sphere \"{} {} {}\" radius 0.35 resolution 12\n",
                                10*p1.x(),10*p1.y(),10*p1.z()
                                );
            // Draw edge cylinder
            out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.3 resolution 12\n",
                       10*p1.x(),10*p1.y(),10*p1.z(),
                       10*p2.x(),10*p2.y(),10*p2.z()
                       );
            /*
            // Draw thin line from central atom to vertexes
            p2 = lip.smoothed_surf_marker;
            out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.1 resolution 12\n",
                                10*p1.x(),10*p1.y(),10*p1.z(),
                                10*p2.x(),10*p2.y(),10*p2.z()
                                );
            */
        }

        // Patch normals
        p1 = lip.surf_marker;
        p2 = p1 + lip.patch.normal*0.5;
        p3 = p1 + lip.patch.normal*0.7;
        out1 += fmt::format("draw color white\n");
        out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2 resolution 12\n",
                   10*p1.x(),10*p1.y(),10*p1.z(),
                   10*p2.x(),10*p2.y(),10*p2.z() );
        out1 += fmt::format("draw cone \"{} {} {}\" \"{} {} {}\" radius 0.3 resolution 12\n",
                   10*p2.x(),10*p2.y(),10*p2.z(),
                   10*p3.x(),10*p3.y(),10*p3.z() );
        // Actual atom
        out1 += fmt::format("draw sphere \"{} {} {}\" radius 0.8 resolution 12\n",
                            10*p1.x(),10*p1.y(),10*p1.z() );

        // fitted normals
        p1 = lip.smoothed_surf_marker;
        out1 += fmt::format("draw color cyan\n");
        p2 = p1 + lip.normal*0.5;
        p3 = p1 + lip.normal*0.7;
        out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.5 resolution 12\n",
                   10*p1.x(),10*p1.y(),10*p1.z(),
                   10*p2.x(),10*p2.y(),10*p2.z() );
        out1 += fmt::format("draw cone \"{} {} {}\" \"{} {} {}\" radius 0.7 resolution 12\n",
                   10*p2.x(),10*p2.y(),10*p2.z(),
                   10*p3.x(),10*p3.y(),10*p3.z() );

        // Smoothed atom
        out1 += fmt::format("draw sphere \"{} {} {}\" radius 1.3 resolution 12\n",
                   10*p1.x(),10*p1.y(),10*p1.z() );


        /*
        // Tails marker
        Vector3f pt = lip.tail_marker;
        out1 += fmt::format("draw color red\n");
        out1 += fmt::format("draw sphere \"{} {} {}\" radius 1.0 resolution 12\n",
                   10*pt.x(),10*pt.y(),10*pt.z() );
        */


        // Inclusion atoms (if any)
        for(size_t i=0; i<lip.surf.active_inclusion_indexes.size(); ++i){
            int ind = lip.surf.active_inclusion_indexes[i];
            // Actual position
            p1 = inclusion.xyz(lip.inclusion_neib[ind]);
            out1 += fmt::format("draw color green\n");
            out1 += fmt::format("draw sphere \"{} {} {}\" radius 1.2 resolution 12\n",
                       10*p1.x(),10*p1.y(),10*p1.z() );
            // Fitted position
            p2 = lip.surf.inclusion_coord.col(ind) + lip.smoothed_surf_marker;
            //p2 = box.closest_image(p2,p1);
            out1 += fmt::format("draw color iceblue\n");
            out1 += fmt::format("draw sphere \"{} {} {}\" radius 1.0 resolution 12\n",
                                10*p2.x(),10*p2.y(),10*p2.z() );
            // Connection line to involved lipid
            p3 = lip.smoothed_surf_marker;
            out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2 resolution 12\n",
                                10*p3.x(),10*p3.y(),10*p3.z(),
                                10*p2.x(),10*p2.y(),10*p2.z() );
        }

        // Principal directions
        // fitted normals
        p1 = lip.smoothed_surf_marker;
        Vector3f dir = lip.surf.principal_directions.col(0);
        p2 = p1 + dir * 0.3; //lip.surf.principal_curvatures(0);
        p3 = p2 + dir * 0.2;
        out1 += fmt::format("draw color blue\n");
        out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.3 resolution 12\n",
                            10*p1.x(),10*p1.y(),10*p1.z(),
                            10*p2.x(),10*p2.y(),10*p2.z() );
        out1 += fmt::format("draw cone \"{} {} {}\" \"{} {} {}\" radius 0.4 resolution 12\n",
                            10*p2.x(),10*p2.y(),10*p2.z(),
                            10*p3.x(),10*p3.y(),10*p3.z() );

        dir = lip.surf.principal_directions.col(1);
        p2 = p1 + dir * 0.3; //lip.surf.principal_curvatures(1);
        p3 = p2 + dir * 0.2;
        out1 += fmt::format("draw color red\n");
        out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.3 resolution 12\n",
                            10*p1.x(),10*p1.y(),10*p1.z(),
                            10*p2.x(),10*p2.y(),10*p2.z() );
        out1 += fmt::format("draw cone \"{} {} {}\" \"{} {} {}\" radius 0.4 resolution 12\n",
                            10*p2.x(),10*p2.y(),10*p2.z(),
                            10*p3.x(),10*p3.y(),10*p3.z() );
    }

    // Output area and normals plots
    make_dir_if_needed(path);
    auto out_all = fmt::output_file(fmt::format("{}/areas_all.tcl",path));
    out_all.print("{}",out1);
    out_all.close();

    for(size_t i=0; i<lipids.size(); ++i){
        all_surf_sel.beta(i) = 10*lipids[i].mean_curvature;
    }
    all_surf_sel.write(fmt::format("{}/areas_all.pdb",path));
}


void LipidMembrane::compute_averages()
{
    // Run post-processing for all groups
    for(auto& gr: groups){
        gr.post_process();
    }
}

void LipidMembrane::write_averages(const string& out_dir)
{        
    string s;
    s += "Run summary:\n";
    s += fmt::format("Lipid species ({}): \n",species.size());
    for(auto& gr: groups){
        s += gr.summary();
    }

    // Print summary
    fmt::print("{}",s);

    // Get output dir path and create output dirs if needed
    auto path = make_dir_if_needed(out_dir);
    log->debug("Output path: {}",path);

    // Save summary to file
    auto out = fmt::output_file((path / "summary.dat").native());
    out.print("{}",s);
    out.close();

    // Write files for species properties
    for(auto const& g: groups){
        g.save_properties_table(path);
        g.save_per_species_properties(path);
    }
}



void pteros::mean_std_from_accumulated(Vector2f &storage, float N){
    if(N>0){
        float s1 = storage[0];
        float s2 = storage[1];
        storage[0] = s1/N;
        storage[1] = sqrt(s2/N - s1*s1/N/N);
    } else {
        storage.fill(0.0);
    }
}
