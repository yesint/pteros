#include "pteros/python/compiled_plugin.h"

class center: public Compiled_plugin_base {
public:
    center(Trajectory_processor* pr, Options_tree* opt): Compiled_plugin_base(pr,opt) {
    }
protected:
    void pre_process(){
        cout << "In compiled_center pre_process" << endl;
    }

    bool process_frame(const Frame_info &info){
        cout << "In compiled_center process_frame " << info.valid_frame << endl;
    }

    void post_process(const Frame_info &info){
        cout << "In compiled_center post_process" << endl;
    }
};

CREATE_COMPILED_PLUGIN(center)
