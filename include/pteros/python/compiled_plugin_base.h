#ifndef COMPILED_PLUGIN_BASE_H
#define COMPILED_PLUGIN_BASE_H

#include "pteros/analysis/consumer.h"
#include "pteros/analysis/trajectory_processor.h"

namespace pteros {

struct Compiled_plugin_base: public Consumer {

    Compiled_plugin_base(Trajectory_processor* proc, Options_tree* opt): Consumer(proc){
        options = opt;
    }

protected:
    /// Options for this particular plugin instance
    Options_tree* options;

    /// Unique prefix for output files for this plugin issue
    /// Set by driver program
    std::string prefix;
};

}

#endif
