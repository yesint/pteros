#pragma once

#include <vector>
#include <Eigen/Core>

namespace pteros {

struct LocalPatch {
    // Id's and distances correspond to each other
    std::vector<int> neib_id;

    // Normal
    Eigen::Vector3f normal;

    // Deduplicate neib_id
    void sort_and_remove_duplicates();
};


}
