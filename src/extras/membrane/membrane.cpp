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



void fit_quad_surface(const MatrixXf& coord,
                      Matrix<float,6,1>& res,
                      Vector3f& smoothed_first_point,
                      float& rms){
    int N = coord.cols();
    // We fit with polynomial fit = A*x^2 + B*y^2 + C*xy + D*x + E*y + F
    // Thus we need a linear system of size 6
    Matrix<float,6,6> m;
    Matrix<float,6,1> rhs; // Right hand side and result

    // Functions returning powers after coefficients (x^2, y^2, x*y, x, y, 1)
    vector<std::function<float(int)> > coef(6);
    coef[0] = [&coord](int j){return coord.col(j)(0)*coord.col(j)(0);};
    coef[1] = [&coord](int j){return coord.col(j)(1)*coord.col(j)(1);};
    coef[2] = [&coord](int j){return coord.col(j)(0)*coord.col(j)(1);};
    coef[3] = [&coord](int j){return coord.col(j)(0);};
    coef[4] = [&coord](int j){return coord.col(j)(1);};
    coef[5] = [&coord](int j){return 1.0;};

    // Now form the matrix
    m.fill(0.0);
    rhs.fill(0.0);

    for(int r=0;r<6;++r){ //rows
        for(int c=0;c<6;++c){ //columns
            m(r,c) = 0.0;
            for(int j=0;j<N;++j){
                m(r,c) += coef[r](j)*coef[c](j);
            }
        }
        // Now rhs
        for(int j=0;j<N;++j){
            rhs(r) += coef[r](j)*coord.col(j)(2);
        }
    }

    // Now solve
    res = m.colPivHouseholderQr().solve(rhs);

    // Compute RMS of fitting
    rms = 0.0;
    for(int j=0;j<N;++j){
        float fit = 0.0;
        for(int r=0;r<6;++r) fit += coef[r](j)*res[r];
        rms += pow(coord.col(j)(2)-fit,2);
    }
    rms = sqrt(rms/float(N));

    // Get smoothed surface point
    smoothed_first_point.fill(0.0);
    for(int r=0;r<6;++r) smoothed_first_point(2) += coef[r](0)*res[r];
}

float compute_area(const Eigen::MatrixXf& coord, double dist, std::vector<int>& neib){
    // Perform Voronoi computation in the tangent plane
    using namespace voro;

    voronoicell_neighbor c;

    // Z dimension of cotainer is 1, so volume=area
    container con(-dist,dist,
                  -dist,dist,
                  -0.5,0.5,
                  10,10,10,false,false,false,8);


    particle_order po;
    // Cycle over particles, first particle is our target
    for(int i=0;i<coord.cols();++i){
        con.put(po, i, coord.col(i)(0), coord.col(i)(1), 0.0); // Z is zero
    }
    c_loop_order clo(con,po);
    if(! clo.start()) return -1.0; // Get first point, which is our target
    if(! con.compute_cell(c,clo) ) return -1.0;

    // Check if any of vertices is on outer bounding box and skip cell if they are
    vector<double> vert;
    c.vertices(0,0,0,vert);
    for(int i=0; i<vert.size(); i+=3){
        //     X                         Y
        if(abs(abs(vert[i])-dist)<=1.0e-5 || abs(abs(vert[i+1])-dist)<=1.0e-5){
            //cout << "| " << vert[i] << " " << vert[i+1] << " " << vert[i+2] << " " << abs(vert[i])-dist << endl;
            return -1.0;
        }
    }

    // Extract neighbors of this particle
    std::vector<int> n;
    c.neighbors(n);
    // Filter out negatives
    neib.clear();
    for(int i: n)
        if(i>0) neib.push_back(i);

/*
    cout << "Npoints: " << coord.cols() << endl;
    // Output the particle positions in gnuplot format
    con.draw_particles("points.gnu");

    // Output the Voronoi cells in gnuplot format
    con.draw_cells_gnuplot("points_v.gnu");
*/
    return c.volume();
}


void get_curvature(const Matrix<float,6,1>& coeffs, float& gaussian_curvature, float& mean_curvature){
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

    float E_ = 1.0+coeffs[3]*coeffs[3];
    float F_ = coeffs[3]*coeffs[4];
    float G_ = 1.0+coeffs[4]*coeffs[4];

    float L_ = 2.0*coeffs[0];
    float M_ = coeffs[2];
    float N_ = 2.0*coeffs[1];

    //Curvatures:
    gaussian_curvature = (L_*N_-M_*M_)/(E_*G_-F_*F_);
    mean_curvature = 0.5*(E_*N_-2.0*F_*M_+G_*L_)/(E_*G_-F_*F_);
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
    neib.clear();
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
    head_marker = head_marker_sel.center(true);
    tail_marker = tail_marker_sel.center(true);
    mid_marker = mid_marker_sel.center(true);

    pos_saved = mid_marker_sel.xyz(0);
    mid_marker_sel.xyz(0) = mid_marker;
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
    for(int i: lip.neib){
        //cout << i << " " << lip.membr_ptr->lipids.size() << endl;
        around[lip.membr_ptr->lipids[i].name] += 1;
    }

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

void LipidMembrane::compute_properties(float d, bool use_external_normal, Vector3f_const_ref external_pivot, Vector3i_const_ref external_dist_dim)
{
    Eigen::Matrix3f tr, tr_inv;

    // Set markers for all lipids
    // This unwraps each lipid
    for(auto& l: lipids) l.set_markers();    

    // Get connectivity
    vector<Vector2i> bon;
    vector<float> dist;
    search_contacts(d,all_mid_sel,bon,dist,false,fullPBC);

    // Convert the list of bonds to convenient form
    // atom ==> 1 2 3...
    vector<vector<int> > conn(all_mid_sel.size());
    for(int i=0;i<bon.size();++i){
        conn[bon[i](0)].push_back(bon[i](1));
        conn[bon[i](1)].push_back(bon[i](0));
    }    

    // Process lipids
    for(int i=0;i<lipids.size();++i){
        LipidMolecule& lip = lipids[i];

        // conn[i] contains local indexes in all_mid_sel
        // Sort it to have predictable order
        std::sort(conn[i].begin(),conn[i].end());

        // Save local selection for this lipid
        lip.local_sel = all_mid_sel.select(conn[i]);

        if(lip.local_sel.size()==0){
            log->warn("Empty locality of lipid {}!",i);
            continue;
        }

        // Local coordinate axes
        Matrix3f axes;

        //-----------------------------------
        // Compute normal
        //-----------------------------------

        Vector3f normal;

        if(use_external_normal){
            // Compute external normal as a vector from pivot to COM of mid_sel (set as marker currently)
            // over specified dimensions
            axes.col(2) = (lip.mid_marker-external_pivot).array()*external_dist_dim.cast<float>().array();

            // Find two vectors perpendicular to normal
            if( axes.col(2).dot(Vector3f(1,0,0)) != 1.0 ){
                axes.col(0) = axes.col(2).cross(Vector3f(1,0,0));
            } else {
                axes.col(0) = axes.col(2).cross(Vector3f(0,1,0));
            }
            axes.col(1) = axes.col(2).cross(axes.col(0));
        } else {
            // Compute normals from inertia axes
            // Create selection for locality of this lipid including itself
            Selection local_self(lip.local_sel);
            local_self.append(all_mid_sel.index(i)); // Add central atom

            // Get inertial axes
            Vector3f moments;
            local_self.inertia(moments,axes,fullPBC); // Have to use periodic variant
            // axes.col(2) will be a normal
        }

        normal = axes.col(2);

        // transformation matrix to local basis
        for(int j=0;j<3;++j) tr.col(j) = axes.col(j).normalized();
        tr_inv = tr.inverse();

        // Need to check direction of the normal
        float ang = angle_between_vectors(normal, lip.head_marker-lip.tail_marker);
        if(ang < M_PI_2){
            lip.normal = normal;
            lip.tilt = rad_to_deg(ang);
        } else {
            lip.normal = -normal;
            lip.tilt = rad_to_deg(M_PI-ang);
        }
        lip.normal.normalized();

        // Create array of local points in local basis
        MatrixXf coord(3,lip.local_sel.size()+1);
        Vector3f c0 = lip.mid_marker; // Real coord of central point - the marker
        coord.col(0) = Vector3f::Zero(); // Local coord of central point is zero
        for(int j=0; j<lip.local_sel.size(); ++j)
            coord.col(j+1) = tr_inv * system->box(0).shortest_vector(c0,lip.local_sel.xyz(j));

        //-----------------------------------
        // Smooth and find local curvatures
        //-----------------------------------

        //if(do_curvature){
            // Fit a quad surface
            Matrix<float,6,1> res;
            Vector3f sm; // smoothed coords of central point
            fit_quad_surface(coord,res,sm,lip.quad_fit_rms);
            // Compute curvatures using quad fit coeefs
            get_curvature(res, lip.gaussian_curvature, lip.mean_curvature);

            // Get smoothed surface point in lab coords
            lip.smoothed_mid_xyz = tr*sm + lip.mid_marker;
        //}

        //-----------------------------------
        // Area and neighbours
        //-----------------------------------
        vector<int> neib;        
        lip.area = compute_area(coord,10.0,neib);
        if(lip.area>0){
            // Save coordination number of the lipid
            lip.coord_number = neib.size();
            // Neighbours are in terms of point indexes in order of their addition
            // this is the same as order in conn[i]. 0 is central atom, neighbors start from 1.
            // conn[i] contains local indexes in all_mid_sel
            lip.neib.resize(neib.size());
            for(int j=0;j<neib.size();++j){
                //cout << j << " " << neib[j] << " " << conn[i].size() <<  " " << conn[i][neib[j]] << endl;
                lip.neib[j] = conn[i][neib[j]-1];
            }
        } else {
            // Something is wrong
            lip.coord_number = -1;
            lip.neib.clear();
        }


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

