#ifndef COMPILED_PLUGIN_BASE_H
#define COMPILED_PLUGIN_BASE_H

#include "pteros/analysis/consumer.h"
#include "pteros/analysis/trajectory_processor.h"

namespace pteros {

struct Compiled_plugin_base: public Consumer {

    Compiled_plugin_base(Trajectory_processor* proc, Options_tree* opt): Consumer(proc){
        options = opt;
    }

    /// Unique textual label for this plugin instance
    /// Set by driver program
    std::string label;

protected:
    /// Options for this particular plugin instance
    Options_tree* options;
};

}

#endif
