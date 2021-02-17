#pragma once

#include "pteros/core/system.h"
#include "pteros/analysis/frame_info.h"

namespace pteros {

/// Auxiliary container for sending data to tasks
struct Data_container {
    /// Frame itself
    Frame frame;
    /// Frame information
    Frame_info frame_info;
};

}
