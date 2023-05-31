#pragma once

#include <vector>
#include <Eigen/Core>

namespace pteros {

struct LocalPatch {
    // Id's and distances correspond to each other
    std::vector<int> neib_id;
    std::vector<float> neib_dist;
    // Local coordinate axes in lab space
    Eigen::Matrix3f axes;
    // Transformations to and from local cordinates
    Eigen::Matrix3f to_lab,to_local;
    // Normal
    Eigen::Vector3f normal;
    // Original central point
    Eigen::Vector3f original_center;

    // Deduplicate neib_id
    void sort_and_remove_duplicates();
};


}
