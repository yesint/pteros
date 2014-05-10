
#include "pteros/pteros.h"
#include "pteros/python/compiled_plugin.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

class PLUGIN_NAME: public Compiled_plugin_base {
public:

    PLUGIN_NAME(Trajectory_processor* pr, const Options& opt): Compiled_plugin_base(pr,opt) { }

    string help(){
        return  "Purpose:\n"
                "\tPut purpose of your plugin here\n"
                "Output:\n"
                "\tDescription of its output\n"
                "Options:\n"
                "\tAny options";
    }

protected:
    void pre_process(){
        string sel_text = options("selection").as_string();
        use_mass = options("mass_weighted","false").as_bool();
        sel.modify(system,sel_text);
        cout << "Working on selection " << sel.get_text() << endl;
        cout << "There are " << sel.size() << " atoms in selection" << endl;
    }

    void process_frame(const Frame_info &info){
        cout << "Frame " << info.absolute_frame << " time " << info.absolute_time << endl;
        cout << "Selection center: " << sel.center(use_mass).transpose() << endl;
    }

    void post_process(const Frame_info &info){
        cout << "Finished" << endl;
    }
private:
    Selection sel;
    bool use_mass;
};

CREATE_COMPILED_PLUGIN(PLUGIN_NAME)
