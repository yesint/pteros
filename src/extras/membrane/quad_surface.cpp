#include "pteros/extras/membrane/quad_surface.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
//#include "voro++.hh"
#include "voro++_2d.hh"
#include <vector>

#include <fmt/ranges.h>

#include "pteros/extras/membrane/lipid_membrane.h"

using namespace pteros;
using namespace std;
using namespace Eigen;

/*
void QuadSurface::fit_quad(){
    int N = lip_coord.cols();
    // We fit with polynomial fit = A*x^2 + B*y^2 + C*xy + D*x + E*y + F
    // Thus we need a linear system of size 6
    Matrix<float,6,6> m;
    Matrix<float,6,1> rhs; // Right hand side and result
    m.fill(0.0);
    rhs.fill(0.0);

    Matrix<float,6,1> powers;
    powers(5) = 1.0; //free term, the same everywhere
    for(int j=0;j<N;++j){ // over points
        powers(0) = lip_coord(0,j)*lip_coord(0,j); //xx
        powers(1) = lip_coord(1,j)*lip_coord(1,j); //yy
        powers(2) = lip_coord(0,j)*lip_coord(1,j); //xy
        powers(3) = lip_coord(0,j); //x
        powers(4) = lip_coord(1,j); //y
        // Fill the matrix
        m.noalias() += powers * powers.transpose();
        // rhs
        rhs += powers*lip_coord(2,j);
    }

    // Now solve
    quad_coefs = m.ldlt().solve(rhs);

    // Set curvatures and normal
    compute_curvature_and_normal();

    // Fit central point
    fitted_central_point(0) = lip_coord(0,0);
    fitted_central_point(1) = lip_coord(1,0);
    fitted_central_point(2) = evalZ(lip_coord(0,0),lip_coord(1,0));

    //fmt::print(">>{} {}\n",fitted_central_point(2),quad_coefs.transpose());
}
*/

void QuadSurface::fit_quad(Vector3f_const_ref normal){
    // Compute transforms
    create_transforms_from_normal(normal);

    // We fit with polynomial fit = A*x^2 + B*y^2 + C*xy + D*x + E*y + F
    // Thus we need a linear system of size 6
    Matrix<float,6,6> m;
    Matrix<float,6,1> rhs; // Right hand side and result
    m.fill(0.0);
    rhs.fill(0.0);

    Matrix<float,6,1> powers;
    powers(5) = 1.0; //free term, the same everywhere
    for(int j=0;j<lip_coord.cols();++j){ // over points
        // Bring point to local basis
        Vector3f loc = to_local*lip_coord.col(j);
        // Compute powers
        powers(0) = loc(0)*loc(0); //xx
        powers(1) = loc(1)*loc(1); //yy
        powers(2) = loc(0)*loc(1); //xy
        powers(3) = loc(0); //x
        powers(4) = loc(1); //y
        // Fill the matrix
        m.noalias() += powers * powers.transpose();
        // rhs
        rhs += powers*loc(2);
    }

    // Now solve
    quad_coefs = m.ldlt().solve(rhs);

    // Set curvatures and normal
    compute_curvature_and_normal();

    // Fit central point
    fitted_central_point = to_lab*Vector3f{0.0,0.0,evalZ(0.0,0.0)};
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
    using namespace voro2d;

    voronoicell_neighbor_2d cell;
    // Initialize with large volume by default in XY and size 1 in Z
    const float side_wall = 100;
    cell.init(-side_wall,side_wall,-side_wall,side_wall);
    // Cut by planes for all points
    for(int i=1; i<lip_coord.cols(); ++i){ // From 1, central point is 0 and ignored
        // Bring point to local basis
        Vector3f loc = to_local*lip_coord.col(i);
        bool ret = cell.nplane(loc(0),loc(1),i); // Pass i as neigbour pid
        if(!ret){ // Fail if the cell is invalid
            in_plane_area = 0;
            surf_area = 0;
            neib_id.clear();
            return;
        }
    }

    // If inclusion is involved add planes from its atoms
    active_inclusion_indexes.clear();
    for(int i=0; i<inclusion_coord.cols(); ++i){
        Vector3f loc = to_local*inclusion_coord.col(i);
        if(abs(loc(2))<inclusion_h_cutoff){
            // Pass large neigbour pid to differentiate
            cell.nplane(loc(0),loc(1),i+10000);
        }
    }

    vector<int> neib_list;
    vector<int> face_vert;
    vector<double> vert_coords;

    cell.neighbors(neib_list);
    cell.vertices_ordered(vert_coords);

    area_vertexes.clear();
    area_vertexes.reserve(neib_list.size());

    // Check if all vertices are inside the box
    for(size_t i=0;i<neib_list.size();++i){
        if(neib_list[i]<0){
            // This is invalid lipid
            in_plane_area = 0;
            surf_area = 0;
            neib_id.clear();
            return;
        }

        // Add vertex
        auto const& x = vert_coords[2*i];
        auto const& y = vert_coords[2*i+1];
        area_vertexes.push_back(Vector3f(x,y,0.0));
    }

    // In-plane area.
    in_plane_area = cell.area();

    // Compute area on surface
    // Project all vertexes to the surface and bring it to lab coordinates
    surf_area = 0.0;
    for(Vector3f& v: area_vertexes) v = to_lab*Vector3f{v(0),v(1),evalZ(v(0),v(1))};

    // Sum up areas of triangles. Points are not duplicated.
    for(int i=0; i<area_vertexes.size(); ++i){
        int ii = (i<area_vertexes.size()-1) ? i+1 : 0;
        surf_area += 0.5*(area_vertexes[i]-fitted_central_point)
                .cross(area_vertexes[ii]-fitted_central_point)
                .norm();
    }

    // Set list of neigbours
    neib_id.clear();
    for(int id: neib_list){
        // Ignore inclusions
        if(id<10000){
            neib_id.push_back(id);
        } else {
            // This is an inclusion
            // Save index of this inclusion atom and project to surface
            int ind = id-10000;
            active_inclusion_indexes.push_back(ind);
            Vector3f v = to_local*inclusion_coord.col(ind);
            inclusion_coord.col(ind) = to_lab*Vector3f{v(0),v(1),evalZ(v(0),v(1))};
        }
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
    fitted_normal = to_lab*Vector3f(D(),E(),-1.0).normalized();
    // Orientation of the normal could be wrong!
    // Have to be flipped according to lipid orientation later

    /* Principal curvatures
        The principal curvatures k1 and k2 are the eigenvalues
        and the principal directions are eigenvectors
        of the shape operator W:
        W = [I]^-1 * [II]
        W = 1/(EG - F^2) * [E L - F M, E M - F N]
                            [G M - F L, G N - F M]
    */
    Matrix2f W;
    W(0,0) = E_*L_ - F_*M_;  W(0,1) = E_*M_ - F_*N_;
    W(1,0) = G_*M_ - F_*L_;  W(1,1) = G_*N_ - F_*M_;
    W *= 1.0/(E_*G_-F_*F_);
    // W is symmetric despite the equations seems to be not!
    Eigen::SelfAdjointEigenSolver<Matrix2f> solver(W);
    principal_directions.col(0) = to_lab*Vector3f{solver.eigenvectors()(0,0),solver.eigenvectors()(1,0),0.0};
    principal_directions.col(1) = to_lab*Vector3f{solver.eigenvectors()(0,1),solver.eigenvectors()(1,1),0.0};
    principal_curvatures = solver.eigenvalues();
}

void QuadSurface::create_transforms_from_normal(Vector3f_const_ref normal)
{
    Matrix3f axes;
    // Set axes as two perpendicular vectors
    axes.col(2) = normal;
    axes.col(0) = axes.col(2).cross(Vector3f{1,0,0});
    axes.col(1) = axes.col(2).cross(axes.col(0));

    // Transforms
    to_lab = axes.colwise().normalized();
    to_local = to_lab.inverse();
}
