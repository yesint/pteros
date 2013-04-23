#include "pteros/python/compiled_plugin.h"

using namespace pteros;

class center: public Compiled_plugin_base {
public:
    center(Trajectory_processor* pr, Options_tree* opt): Compiled_plugin_base(pr,opt) {
    }
protected:
    void pre_process(){
        string sel_text = options->get_value<string>("selection");
        use_mass = options->get_value<bool>("mass_weighted",false);
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

CREATE_COMPILED_PLUGIN(center)
