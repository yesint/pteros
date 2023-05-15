#pragma once

#include <Eigen/Geometry>
#include "pteros/core/selection.h"

namespace pteros {

Eigen::Affine3f fit_transform(const Selection& sel1, const Selection& sel2, bool translate_to_zero=true);
Eigen::MatrixXf rmsd_matrix(Selection& sel);

}
