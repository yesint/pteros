#pragma once

namespace pteros {

// Properties of individual lipids
struct LipidProperties {
    Eigen::Vector3f normal;
    float tilt;
    float area;
    int coord_number;
    float gaussian_curvature;
    float mean_curvature;
};

}
