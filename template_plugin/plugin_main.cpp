
#include "pteros/pteros.h"
#include "pteros/python/compiled_plugin.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

TASK_SERIAL(PLUGIN_NAME)
public:

    string help() override {
        return  "Purpose:\n"
                "\tPut purpose of your plugin here\n"
                "Output:\n"
                "\tDescription of its output\n"
                "Options:\n"
                "\tAny options";
    }

protected:
    void pre_process() override {
        string sel_text = options("selection").as_string();
        use_mass = options("mass_weighted","false").as_bool();
        sel.modify(system,sel_text);
        log->info("Working on selection '{}'", sel.get_text());
        log->info("There are {} atoms in selection", sel.size());
    }

    void process_frame(const FrameInfo &info) override {
        log->info("Frame: {}, time: {}, center: {}", info.absolute_frame, info.absolute_time, sel.center(use_mass).transpose());
    }

    void post_process(const FrameInfo &info) override {
        log->info("Finished");
    }
private:
    Selection sel;
    bool use_mass;
};

CREATE_COMPILED_PLUGIN(PLUGIN_NAME)
