#include "pteros/extras/membrane/quad_surface.h"
#include <Eigen/Dense>
#include "voro++.hh"
#include <vector>

#include <fmt/ranges.h>

using namespace pteros;
using namespace std;
using namespace Eigen;

void QuadSurface::fit_to_points(MatrixXf_const_ref coord){
    int N = coord.cols();
    // We fit with polynomial fit = A*x^2 + B*y^2 + C*xy + D*x + E*y + F
    // Thus we need a linear system of size 6
    Matrix<float,6,6> m;
    Matrix<float,6,1> rhs; // Right hand side and result
    m.fill(0.0);
    rhs.fill(0.0);

    Matrix<float,6,1> powers;
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

void QuadSurface::project_point_to_surface(Vector3f_ref p){
    p(2) =  evalZ(p(0),p(1));
}

float QuadSurface::evalZ(float x, float y){
    return    A()*x*x
            + B()*y*y
            + C()*x*y
            + D()*x
            + E()*y
            + F();
}

void QuadSurface::compute_voronoi(float inclusion_h_cutoff){
    // The center of the cell is assumed to be at {0,0,0}
    // First point is assumed to be the central one and thus not used
    using namespace voro;

    voronoicell_neighbor cell;
    // Initialize with large volume by default in XY and size 1 in Z
    const float side_wall = 100;
    cell.init(-side_wall,side_wall,-side_wall,side_wall,-0.5,0.5);
    // Cut by planes for all points
    for(int i=1; i<fitted_points.cols(); ++i){ // From 1, central point is 0 and ignored
        cell.nplane(fitted_points(0,i),fitted_points(1,i),0.0,i); // Pass i as neigbour pid
    }

    // If inclusion is involved add planes from its atoms
    for(int i=0; i<inclusion_coord.cols(); ++i){
        if(abs(inclusion_coord(2,i))<inclusion_h_cutoff){
            cell.nplane(inclusion_coord(0,i),inclusion_coord(1,i),0.0,i*10000);
            // Pass large neigbour pid to differenciate
        }
    }

    vector<int> neib_list;
    vector<int> face_vert;
    vector<double> vert_coords;

    //cell.draw_gnuplot(0,0,0,"cell.dat");

    cell.neighbors(neib_list);
    //fmt::print("{}\n",fmt::join(neib_list, ", "));

    // The top and bottom walls have indexes -5 and -6
    // Side walls has indexes -1..-4

    cell.face_vertices(face_vert);
    // In face_vert the first entry is a number k corresponding to the number
    // of vertices making up a face,
    // this is followed by k additional entries describing which vertices make up this face.
    // For example, (3, 16, 20, 13) is a triangular face linking vertices 16, 20, and 13

    cell.vertices(0,0,0,vert_coords);
    //vert_coords this corresponds to a vector of triplets (x, y, z)
    // describing the position of each vertex.

    // Collect verteces which make up cell area in XY plane
    int vert_ind;
    float x,y;
    area_vertexes.clear();
    area_vertexes.reserve(cell.number_of_faces()-2); // 2 faces are top and bottom walls

    // We need to take only face at the upped wall
    bool invalid_lipid = false;
    int j = 0;
    for(int i=0;i<neib_list.size();++i){
        if(neib_list[i]==-5){ // Only take the face at upper wall (-5)
            // Cycle over vertices of this face
            for(int ind=0; ind<face_vert[j]; ++ind){
                vert_ind = 3*face_vert[j+1+ind];
                x = vert_coords[vert_ind];
                y = vert_coords[vert_ind+1];
                // We are not interested in Z component
                if(abs(x)<0.9*side_wall && abs(y)<0.9*side_wall){
                    // Add vertex
                    area_vertexes.push_back(Vector3f(x,y,0.0));
                } else {
                    // This is invalid lipid
                    invalid_lipid = true;
                    break; // Bail out
                }
            }
            // We are done
            break;
        }
        j += face_vert[j]+1; // Next entry in face_vert array
    }

    // In case of invalid lipid give up and return early
    if(invalid_lipid){
        in_plane_area = 0;
        surf_area = 0;
        neib_id.clear();
        return;
    }

    // In-plane area. Since Z size of the cell is 1 volume=area
    in_plane_area = cell.volume();

    // Compute area on surface
    // Project all vertexes to the surface
    surf_area = 0.0;
    for(Vector3f& v: area_vertexes) project_point_to_surface(v);
    // Sum up areas of triangles. Points are not duplicated.
    for(int i=0; i<area_vertexes.size(); ++i){
        int ii = (i<area_vertexes.size()-1) ? i+1 : 0;
        surf_area += 0.5*(area_vertexes[i]-fitted_points.col(0))
                .cross(area_vertexes[ii]-fitted_points.col(0))
                .norm();
    }

    // Set list of neigbours
    // Filter out negatives
    neib_id.clear();
    for(int id: neib_list){
        // Ignore walls and inclusions
        if(id>=0 && id<10000) neib_id.push_back(id);
    }
}


void QuadSurface::compute_curvature_and_normal(){
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
