#pragma once

#include <Eigen/Core>
#include "pteros/core/typedefs.h"

namespace pteros {
// Helper class representing a quadric surface for fitting
// Surface is assumed to be centered around point {0,0,0}
// We fit with polynomial fit = A*x^2 + B*y^2 + C*xy + D*x + E*y + F
class QuadSurface {
public:
    // Compute fitting surface given the points in local basis
    void fit_quad();

    // Projects point from the XY plane to the surface
    // existing Z coordinate will be substituted in place by the surface Z value
    inline void project_point_to_surface(Vector3f_ref p);

    // Compute Z point of the surface
    inline float evalZ(float x, float y);

    // Computes voronoi tesselation in the local tangent plane
    // Sets vertexes of points in plane making up the area
    // Computes in-plane area
    void compute_voronoi(float inclusion_h_cutoff);

    void compute_curvature_and_normal();

    inline float A(){ return quad_coefs[0]; }
    inline float B(){ return quad_coefs[1]; }
    inline float C(){ return quad_coefs[2]; }
    inline float D(){ return quad_coefs[3]; }
    inline float E(){ return quad_coefs[4]; }
    inline float F(){ return quad_coefs[5]; }

//private:
    // Coefficients of the quadric surface A,B,C,D,E,F
    Eigen::Matrix<float,6,1> quad_coefs;
    // Fitted normal (unoriented!)
    Eigen::Vector3f fitted_normal;
    // Vertexes for area calculations
    std::vector<Eigen::Vector3f> area_vertexes;
    // List of neighbour ids
    std::vector<int> neib_id;
    // Computed properties
    float fit_rms;
    float in_plane_area;
    float surf_area;
    float gaussian_curvature;
    float mean_curvature;

    // Local coordinates of lipid markers to account for
    Eigen::MatrixXf lip_coord;
    // Local coordinates of inclusion atoms to account for
    Eigen::MatrixXf inclusion_coord;

    Eigen::Vector3f fitted_central_point;
};

}
