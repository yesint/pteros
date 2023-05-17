#pragma once

#include <vector>
#include <Eigen/Core>
#include "pteros/core/selection.h"

namespace pteros {

using namespace Eigen;

class Subset {
public:

    Subset(Selection& sel){
        index = sel.get_index_ptr();
        buffer = &sel.get_system()->frame(sel.get_frame()).coord;
    }

    template <typename D>
    DenseBase<D>& data() {
        return Map<MatrixXf>((float*) buffer->data(),3,index->size())(indexing::all,*index);
    }

    Vector3f& pos(size_t i) {
        return (*buffer)[i];
    }

    size_t size() const {return index->size();}

private:
    std::vector<int>* index;
    std::vector<Vector3f>* buffer;
};

}
