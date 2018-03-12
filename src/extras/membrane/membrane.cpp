/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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



Membrane::Membrane(System *sys, const std::vector<Lipid_descr> &species): system(sys), lipid_species(species)
{
    log = create_logger("membrane");

    Selection all_mid_sel(*system);

    // Creating selections
    for(auto& sp: lipid_species){
        vector<Selection> res;
        system->select(sp.whole_sel_str).split_by_residue(res);
        log->info("Lipid {}: {}", sp.name,res.size());
        for(auto& lip: res){
            auto mol = Lipid(lip,sp);
            mol.set_markers();
            all_mid_sel.append(mol.mid_sel.index(0));
            lipids.push_back(mol);
            index_map[mol.mid_sel.index(0)] = lipids.size()-1;
        }
    }

    // Compute connectivity    
    std::vector<Selection> leafs;
    all_mid_sel.split_by_connectivity(2.0,leafs,true);
    leaflets.resize(leafs.size());
    leaflets_sel.resize(leafs.size());
    // Set lipid leaflets
    for(int i=0;i<leafs.size();++i){
        leaflets_sel[i].set_system(*system);
        for(int a=0; a<leafs[i].size(); ++a){
            int ind = leafs[i].index(a);
            lipids[index_map[ind]].leaflet = i;
            leaflets[i].push_back(index_map[ind]);
            leaflets_sel[i].append(ind);
        }
    }

    // Print statictics
    log->info("Total number of lipids: {}",lipids.size());
    log->info("Number of leaflets: {}",leafs.size());    
}

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

void Membrane::compute_properties(float d, bool use_external_normal, Vector3f_const_ref external_pivot, Vector3i_const_ref external_dist_dim)
{
    // Clear everything
    neighbor_pairs.clear();

    Eigen::Matrix3f tr, tr_inv;

    // Set markers for all lipids
    // This unwraps each lipid
    for(auto& l: lipids) l.set_markers();

    // Compute per leaflet properties
    for(int l=0;l<leaflets.size();++l){
        // Get connectivity in this leaflet
        vector<Vector2i> bon;
        search_contacts(d,leaflets_sel[l],bon,false,true);        

        // Convert the list of bonds to convenient form
        // atom ==> 1 2 3...
        vector<vector<int> > conn(leaflets_sel[l].size());
        for(int i=0;i<bon.size();++i){
            conn[bon[i](0)].push_back(bon[i](1));
            conn[bon[i](1)].push_back(bon[i](0));
        }               

        // Process lipids from this leaflet
        for(int i=0;i<leaflets[l].size();++i){
            // Find current lipid
            Lipid& lip = lipids[leaflets[l][i]];

            // Save local selection for this lipid
            lip.local_sel = leaflets_sel[l](conn[i]);

            if(lip.local_sel.size()==0){
                log->warn("Empty locality of lipid {} in leaflet {}! Skipped.",i,l);
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

                local_self.append(leaflets_sel[l].index(i)); // Add central atom

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
                lip.tilt = ang;
            } else {
                lip.normal = -normal;
                lip.tilt = M_PI-ang;
            }
            lip.normal.normalized();

            //-----------------------------------
            // Smooth and find local curvatures
            //-----------------------------------

            // Create array of local points in local basis
            MatrixXf coord(3,lip.local_sel.size()+1);
            Vector3f c0 = lip.mid_marker; // Coord of central point - the marker
            coord.col(0) = Vector3f::Zero();
            for(int j=0; j<lip.local_sel.size(); ++j)
                coord.col(j+1) = tr_inv * system->box(0).shortest_vector(c0,lip.local_sel.xyz(j));

            // Fit a quad surface
            Matrix<float,6,1> res;
            Vector3f sm; // smoothed coords of central point
            fit_quad_surface(coord,res,sm,lip.quad_fit_rms);
            // Compute curvatures using quad fit coeefs
            get_curvature(res, lip.gaussian_curvature, lip.mean_curvature);

            // Get smoothed surface point in lab coords
            lip.smoothed_mid_xyz = tr*sm + lip.mid_marker;

            //-----------------------------------
            // Area and neighbours
            //-----------------------------------

            vector<int> neib;
            lip.area = compute_area(coord,10.0,neib);
            // Save coordination number of the lipid
            lip.coord_number = neib.size();

            // Add neighbor pairs. Only use nearest neighbor lipids from area computation
            // use only i<j to avoid adding duplicates
            // If area is -1 skip since this lipid has weird surrounding
            if(lip.area>0){
                int cur_ind = index_map[leaflets_sel[l].index(i)];
                for(int j=0; j<neib.size(); ++j){
                    int n = index_map[lip.local_sel.index(neib[j]-1)];
                    if(cur_ind<n) neighbor_pairs.push_back(Vector2i(cur_ind,n));
                }
            }

        } // Over lipids in leaflet

    } // over leaflets

    //-----------------------------------
    // Splay and triangilation
    //-----------------------------------

    // Go over neighbor pairs and compute mean splay
    // Also form neigbhor array as i ==> 1,2,3...
    splay.resize(neighbor_pairs.size());
    neighbors.resize(lipids.size());
    for(int i=0;i<neighbor_pairs.size();++i){
        neighbors[neighbor_pairs[i](0)].push_back(neighbor_pairs[i](1));
        neighbors[neighbor_pairs[i](1)].push_back(neighbor_pairs[i](0));
        // We have to use periodic distance since atoms could be from different images
        Lipid& lip1 = lipids[neighbor_pairs[i](0)];
        Lipid& lip2 = lipids[neighbor_pairs[i](1)];
        Vector3f& n1 = lip1.normal;
        Vector3f& n2 = lip2.normal;

        if(n1.dot(n2)<0) continue;

        auto x= system->box(0).shortest_vector(lip1.mid_marker, lip2.mid_marker);
        float d = x.norm();
        x = x-x.dot(n1)*n1;
        x.normalize();
        Vector3f v1 = (lip1.head_marker-lip1.tail_marker).normalized();
        Vector3f v2 = (lip2.head_marker-lip2.tail_marker).normalized();

        splay[i] = {
                    neighbor_pairs[i](0),
                    neighbor_pairs[i](1),
                    (v2.dot(x)-v1.dot(x)-n2.dot(x)+n1.dot(x))/d
                   };
    }

    //-----------------------------------
    // Order parameter
    //-----------------------------------

    for(auto& lip: lipids){
        // Compute Sz order parameter if the tails are provided
        for(int t=0; t<lip.tail_carbon_indexes.size(); ++t){
            // Go over atoms in tail t
            for(int at=1; at<lip.tail_carbon_indexes[t].size()-1; ++at){
                // Vector from at+1 to at-1
                auto coord1 = system->xyz(lip.tail_carbon_indexes[t][at+1]);
                auto coord2 = system->xyz(lip.tail_carbon_indexes[t][at-1]);
                float ang = angle_between_vectors(coord1-coord2,lip.normal);
                lip.order[t][at-1] = 1.5*pow(cos(ang),2)-0.5;
            }
        }
    }

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

void Membrane::write_vmd_arrows(const string &fname)
{
    ofstream f(fname);
    for(auto& lip: lipids){
        f << tcl_arrow(lip.mid_sel.center(true),lip.mid_sel.center(true)+lip.normal,0.2,"green");
        f << tcl_arrow(lip.tail_sel.center(true),lip.head_sel.center(true),0.2,"red");
    }
    f.close();

}

void Membrane::write_smoothed(const string& fname){
    System out;
    for(auto& l: lipids){
        auto s = out.append(l.mid_sel(0,0));
        s.name(0) = "M";
        s = out.append(l.head_sel(0,0));
        s.xyz(0) = l.smoothed_mid_xyz;
        s.name(0) = "S";
    }
    out().write(fname);
}


//------------------------------------------------------------------------

Lipid::Lipid(const Selection &sel, const Lipid_descr &descr){
    name = descr.name;
    whole_sel = sel;
    head_sel = whole_sel(descr.head_sel_str);
    tail_sel = whole_sel(descr.tail_sel_str);
    mid_sel = whole_sel(descr.mid_sel_str);
    // Fill tail indexes if any
    tail_carbon_indexes.resize(descr.tail_carbon_sels.size());
    order.resize(descr.tail_carbon_sels.size());
    for(int t=0; t<descr.tail_carbon_sels.size(); ++t){
        tail_carbon_indexes[t] = whole_sel(descr.tail_carbon_sels[t]).get_index();
        // Allocate array for order
        order[t].resize(tail_carbon_indexes[t].size()-2);
    }
}

void Lipid::set_markers()
{
    // Unwrap this lipid with leading index of position[0]    
    whole_sel.unwrap(fullPBC, mid_sel.index(0)-whole_sel.index(0));

    // Set markers to COM
    head_marker = head_sel.center(true);
    tail_marker = tail_sel.center(true);
    mid_marker = mid_sel.center(true);
}




