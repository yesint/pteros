#pragma once

#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>

namespace pteros {

class Topmatch {
public:
    Topmatch(){}

    void compute_mapping(const Selection& src, const Selection& target);
    void map_and_transform(Selection& src, Selection& target){}

    virtual ~Topmatch(){}

    std::vector<Eigen::Vector2i> get_mapping(){}
    int map_src_to_target(int i){}
    int map_target_to_src(int i){}

protected:

};

}
