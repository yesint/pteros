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

#include "pteros/extras/membrane/lipid_membrane.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>
#include <fstream>
#include <unordered_set>
#include "voro++.hh"
#include <omp.h>
#include <filesystem>

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

LipidMembrane::LipidMembrane(const System *sys,
                             int ngroups,
                             const vector<LipidSpecies>& sp_list,
                             const Selection &incl,
                             float incl_h_cutoff)
{
    log = create_logger("membrane");
    system_ptr = sys;

    // Add all lipid species provided
    for(auto& sp: sp_list) register_lipid_species(sp);

    // Create groups (they see all registered species)
    groups.reserve(ngroups);
    for(int i=0;i<ngroups;++i) groups.emplace_back(this,i);
    log->info("{} groups created",ngroups);

    // Make local copies of inclusion selection
    inclusion = incl;
    inclusion_h_cutoff = incl_h_cutoff;

    // Print statictics
    log->info("Total number of lipids: {}",lipids.size());

    // Create selection for all surf markesr
    all_surf_sel.set_system(*system_ptr);
    vector<int> ind;
    ind.reserve(lipids.size());
    for(auto& lip: lipids) ind.push_back(lip.surf_marker_sel.index(0));
    all_surf_sel.modify(ind);
}


void LipidMembrane::register_lipid_species(const LipidSpecies &sp)
{
    log->info("Adding lipid species {}",sp.name);
    species.push_back(sp);

    vector<Selection> res;
    system_ptr->select(sp.whole_sel_str).split_by_residue(res);
    log->info("There are {} {} lipids in the system", res.size(),sp.name);
    // Add lipids
    int id = lipids.size();
    for(int i=0; i<res.size(); ++i){
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


void LipidMembrane::compute_properties(float d, float incl_d, OrderType order_type)
{    
    // Set markers for all lipids
    // This unwraps each lipid
    #pragma omp parallel for if (lipids.size() >= 100)
    for(auto& l: lipids) l.set_markers();

    // Get connectivity
    vector<Vector2i> bon;
    vector<float> dist;
    search_contacts(d,all_surf_sel,bon,dist,false,fullPBC);

    // Clear patches for all lipids
    #pragma omp parallel for if (lipids.size() >= 100)
    for(auto& l: lipids) {
        l.patch.neib_id.clear();
        l.patch.neib_dist.clear();
    }

    // Fill patches with id's and distances
    for(size_t i=0;i<bon.size();++i){
        int l1 = bon[i](0);
        int l2 = bon[i](1);
        lipids[l1].patch.neib_id.push_back(l2);
        lipids[l1].patch.neib_dist.push_back(dist[l2]);
        lipids[l2].patch.neib_id.push_back(l1);
        lipids[l2].patch.neib_dist.push_back(dist[l1]);
    }

    // If inclusions are present compute contacts with them
    if(inclusion.size()>0){
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

    #pragma omp parallel for if (lipids.size() >= 100)
    // Compute local coordinates and approximate normals of patches
    for(size_t i=0; i<lipids.size(); ++i){
        auto& patch = lipids[i].patch;
        // Save central point
        patch.original_center = lipids[i].surf_marker;
        // Set local selection for this lipid
        lipids[i].local_sel = all_surf_sel.select(patch.neib_id);
        lipids[i].local_sel_with_self = lipids[i].local_sel;
        lipids[i].local_sel_with_self.append(all_surf_sel.index(i));

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

        // We only include lipids which close enough
        // The cutoff is larger for lipids near inclusions
        float dist_cutoff = lipids[i].inclusion_neib.empty() ? 1.0 : d;
        // Angle tolerance is smaller near inclusions
        float ang_tol = lipids[i].inclusion_neib.empty() ? M_PI_4 : 0.5*M_PI_4;

        for(size_t j=0; j<lipids[i].patch.neib_id.size(); ++j){
            int jl = lipids[i].patch.neib_id[j];
            float d = lipids[i].patch.neib_dist[j];
            const Vector3f& normal2 = lipids[jl].patch.normal;

            if( d<dist_cutoff){
                if(angle_between_vectors(normal1,normal2)>ang_tol) ++n_bad;
                aver_closest += normal2;
            }
        }

        if(n_bad>2){
            // Set to average of closest normals
            log->debug("Trying to fix bad normal for lipid {}", i);
            lipids[i].patch.normal = aver_closest.normalized();
            // Get new tangent axes
            lipids[i].patch.axes.col(0) = lipids[i].patch.normal.cross(lipids[i].patch.axes.col(1));
            lipids[i].patch.axes.col(1) = lipids[i].patch.normal.cross(lipids[i].patch.axes.col(0));
            // Recompute transforms
            lipids[i].patch.to_lab = lipids[i].patch.axes.colwise().normalized();
            lipids[i].patch.to_local = lipids[i].patch.to_lab.inverse();
        }
    }

    //================
    // Process lipids
    //================
    #pragma omp parallel for if (lipids.size() >= 100)
    for(size_t i=0;i<lipids.size();++i){
        LipidMolecule& lip = lipids[i];

        // For lipids with inclusions add neighbors of neigbors
        // in order to increase the patch coverage and compensate
        // for absent partners from inclusion side
        if(!lip.inclusion_neib.empty()){
            auto cur_neibs = lip.patch.neib_id;
            for(int ind: cur_neibs){
                for(size_t i=0; i<lipids[ind].patch.neib_id.size(); ++i){
                    lip.patch.neib_id.push_back(lipids[ind].patch.neib_id[i]);
                    lip.patch.neib_dist.push_back(lipids[ind].patch.neib_dist[i]);
                }
            }
            // Update sellections
            lip.local_sel = all_surf_sel.select(lip.patch.neib_id);
            lip.local_sel_with_self = lip.local_sel;
            lip.local_sel_with_self.append(all_surf_sel.index(i));
        }
        // end inclusions

        // Sort neib_id for correct index match with selection
        sort(lip.patch.neib_id.begin(),lip.patch.neib_id.end());
        // Remove duplicates
        vector<int>::iterator it = unique(lip.patch.neib_id.begin(), lip.patch.neib_id.end());
        lip.patch.neib_id.resize(it - lip.patch.neib_id.begin());

        // Create array of local points in local basis
        MatrixXf coord(3,lip.local_sel.size()+1);
        coord.col(0) = Vector3f::Zero(); // Local coord of central point is zero and it goes first
        for(int j=0; j<lip.local_sel.size(); ++j){
            coord.col(j+1) = lip.patch.to_local
                    * system_ptr->box(0).shortest_vector(lip.surf_marker,
                                                     lip.local_sel.xyz(j));
        }

        // Fit a quadric surface and find local curvatures
        lip.surf.fit_to_points(coord);
        // Set smoothed point
        lip.smoothed_surf_marker_xyz = lip.patch.to_lab*lip.surf.fitted_points.col(0) + lip.surf_marker;

        // If inclusion is present bring neibouring inclusion atoms to local basis
        if(!lip.inclusion_neib.empty()){
            lip.surf.inclusion_coord.resize(3,lip.inclusion_neib.size());
            for(size_t i=0; i<lip.inclusion_neib.size(); ++i){
                lip.surf.inclusion_coord.col(i) = lip.patch.to_local * system_ptr->box(0)
                        .shortest_vector(lip.surf_marker, inclusion.xyz(lip.inclusion_neib[i]));

            }
        }

        lip.surf.compute_voronoi(inclusion_h_cutoff);
        // area
        lip.area = lip.surf.surf_area;

        lip.surf.compute_curvature_and_normal();

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
        // Curvatures
        lip.mean_curvature = lip.surf.mean_curvature;
        lip.gaussian_curvature = lip.surf.gaussian_curvature;
        // Normal
        lip.normal = lip.patch.to_lab*lip.surf.fitted_normal;
        // Tilt
        lip.tilt = rad_to_deg(angle_between_vectors(lip.normal,lip.tail_head_vector));        

        // Set neigbours for lipid
        lip.neib.clear();
        for(size_t i=0; i<lip.surf.neib_id.size(); ++i){
            // Indexes in lip.surf.neib_id start from 1, because 0 is the central point
            // Thus substract 1 and we get a local selection index
            // which are in turn correspond to lipid indexes in all_surf_sel.
            // The neighbours are NOT arranged in triangulated order!
            int ind = lip.surf.neib_id[i]-1;
            lip.neib.push_back( lip.patch.neib_id[ind] );
        }

        // Coordination number
        lip.coord_number = lip.neib.size();

    } // for lipids

    // Unset markers. This restores all correct atomic coordinates for analysis
    for(auto& lip: lipids) lip.unset_markers();

    // Analysis whith requires restored coordinates of markers

    #pragma omp parallel for if (lipids.size() >= 100)
    for(auto& lip: lipids){
        // Tail properties        
        for(auto& t: lip.tails){
            t.compute_order_and_dihedrals(lip.whole_sel,
                                          lip.normal,
                                          order_type);
        }
    }

    // Process groups
    for(auto& gr: groups){
        gr.process_frame();
    }
}

MatrixXf LipidMembrane::get_average_curvatures(int lipid, int n_shells)
{
    MatrixXf m(n_shells,2);
    m.fill(0.0);

    vector<int> neib_n;
    //lipids[i].gaussian_curvature_aver.resize(4);
    //lipids[i].mean_curvature_aver.resize(4);
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
    for(size_t i1=0; i1<lipids.size(); ++i1){
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
                    Vector3f v1 = box.shortest_vector(lipids[i1].smoothed_surf_marker_xyz,lipids[i2].smoothed_surf_marker_xyz);
                    Vector3f v2 = box.shortest_vector(lipids[i1].smoothed_surf_marker_xyz,lipids[i3].smoothed_surf_marker_xyz);
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

    for(int smooth_level=0;smooth_level<5;++smooth_level){
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
            Vector3f p1 = lipids[i].smoothed_surf_marker_xyz;
            s+=fmt::format("draw color {}\n",color_ind[i]);
            s+=fmt::format("draw sphere \"{} {} {}\" radius 1.3 resolution 12\n",
                           p1(0)*10,p1(1)*10,p1(2)*10);

            s+=fmt::format("draw color black\n");
            for(int nb: lipids[i].neib){
                Vector3f p2 = box.closest_image(lipids[nb].smoothed_surf_marker_xyz,p1);
                s+=fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.1\n",
                               p1(0)*10,p1(1)*10,p1(2)*10,
                               p2(0)*10,p2(1)*10,p2(2)*10
                              );
            }

        }



        // Output colored triangles
        for(auto t: triangles){
            Vector3f p1 = lipids[t(0)].smoothed_surf_marker_xyz;
            Vector3f p2 = lipids[t(1)].smoothed_surf_marker_xyz;
            Vector3f p3 = lipids[t(2)].smoothed_surf_marker_xyz;
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

        ofstream f(fmt::format("triangulated_smooth_level_{}.tcl",smooth_level));
        f << s;
        f.close();
    } //smooth level
    //exit(1);
}


void LipidMembrane::write_vmd_visualization(const string &path){
    string out1;
    for(const auto& lip: lipids){
        Vector3f p1,p2,p3;        
        // Area vertices in lab coordinates
        out1+="draw materials on\n";
        out1+="draw material AOEdgy\n";
        out1 += fmt::format("draw color orange\n");
        for(size_t j=0; j<lip.surf.area_vertexes.size(); ++j){
            int j2 = j+1;
            if(j==lip.surf.area_vertexes.size()-1) j2=0;
            p1 = lip.patch.to_lab *lip.surf.area_vertexes[j] + lip.patch.original_center;
            p2 = lip.patch.to_lab *lip.surf.area_vertexes[j2] +lip.patch.original_center;
            out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.3 resolution 12\n",
                       10*p1.x(),10*p1.y(),10*p1.z(),
                       10*p2.x(),10*p2.y(),10*p2.z()
                       );
        }
        // Normal vectors
        //p1 = lip.patch.to_lab*fp + lip.patch.original_center;

        // Patch normals
        p1 = lip.patch.original_center;
        p2 = p1 + lip.patch.normal*1.0;
        p3 = p1 + lip.patch.normal*1.2;
        out1 += fmt::format("draw color white\n");
        out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.2 resolution 12\n",
                   10*p1.x(),10*p1.y(),10*p1.z(),
                   10*p2.x(),10*p2.y(),10*p2.z() );
        out1 += fmt::format("draw cone \"{} {} {}\" \"{} {} {}\" radius 0.3 resolution 12\n",
                   10*p2.x(),10*p2.y(),10*p2.z(),
                   10*p3.x(),10*p3.y(),10*p3.z() );

        // fitted normals
        p1 = lip.smoothed_surf_marker_xyz;
        out1 += fmt::format("draw color cyan\n");
        p2 = p1 + lip.normal*0.75;
        p3 = p1 + lip.normal*1.0;
        out1 += fmt::format("draw cylinder \"{} {} {}\" \"{} {} {}\" radius 0.5 resolution 12\n",
                   10*p1.x(),10*p1.y(),10*p1.z(),
                   10*p2.x(),10*p2.y(),10*p2.z() );
        out1 += fmt::format("draw cone \"{} {} {}\" \"{} {} {}\" radius 0.7 resolution 12\n",
                   10*p2.x(),10*p2.y(),10*p2.z(),
                   10*p3.x(),10*p3.y(),10*p3.z() );

        // Smoothed atom
        out1 += fmt::format("draw sphere \"{} {} {}\" radius 1.5 resolution 12\n",
                   10*p1.x(),10*p1.y(),10*p1.z() );

        // Inclusion atoms (if any)
        out1 += fmt::format("draw color green\n");
        for(size_t i=0; i<lip.inclusion_neib.size(); ++i){
            p1 = inclusion.xyz(lip.inclusion_neib[i]);
            out1 += fmt::format("draw sphere \"{} {} {}\" radius 0.3 resolution 12\n",
                       10*p1.x(),10*p1.y(),10*p1.z() );
        }
    }

    // Output area and normals plots
    ofstream out_all(fmt::format("{}/areas_all.tcl",path));
    out_all << out1;
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

void LipidMembrane::write_averages(string path)
{        
    string s;
    s += "Run summary:\n";
    s += fmt::format("Lipid species ({}): \n",species.size());
    for(auto& gr: groups){
        s += gr.summary();
    }

    // Print summary
    cout << s << endl;

    // Create output path if needed
    if(!std::filesystem::exists(path)){
        log->info("Creating output path {}",path);
        std::error_code ec;
        bool ok = std::filesystem::create_directories(path,ec);
        if(!ok){
            throw PterosError("Unable to create directory {}: {}",path,ec.message());
        }
    }

    // Save summary to file
    ofstream out(path+"/summary.dat");
    out << s;
    out.close();

    // Write files for species properties
    for(size_t g=0;g<groups.size();++g){
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
                // Output order histogram
                //sp.second.order_hist.save_to_file(file_prefix+"order_hist.dat");
                // Around
                sp.second.save_around_to_file(file_prefix+"around.dat");
            }
        }
    }
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
