#pragma once

#include <vector>
#include <Eigen/Core>
#include "pteros/core/selection.h"

namespace pteros {

using namespace Eigen;

/// Thin wrapper, which encapsulates the set of coordinates
/// (or velocities, or forces)
/// and the index of selected positions in this set
class Subset {
public:

    /// fr=-1 means current frame
    Subset(Selection& sel, int fr=-1) {
        index = sel.get_index_ptr();
        size_t wf = fr<0 ? sel.get_frame() : fr;
        buffer = &sel.get_system()->frame(wf).coord;
    }

    template <typename D>
    DenseBase<D>& data() const {
        return Map<MatrixXf>((float*)buffer->data(),3,index->size())(indexing::all,*index);
    }

    Vector3f& pos(size_t i) const {
        return (*buffer)[(*index)[i]];
    }

    size_t size() const {return index->size();}

private:
    std::vector<int>* index;
    std::vector<Vector3f>* buffer;
};

}
