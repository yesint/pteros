#pragma once

#include "pteros/core/system.h"
#include "pteros/analysis/frame_info.h"

namespace pteros {

/// Auxiliary container for sending data to tasks
struct DataContainer {
    /// Frame itself
    Frame frame;
    /// Frame information
    FrameInfo frame_info;
};

}
