#ifndef DATA_CONTAINER_H
#define DATA_CONTAINER_H

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

#endif // DATA_CONTAINER_H
