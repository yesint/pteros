#include "pteros/python/compiled_plugin.h"
#include <fstream>

using namespace pteros;

class center: public Compiled_plugin_base {
public:
    center(Trajectory_processor* pr, const Options& opt): Compiled_plugin_base(pr,opt) {
    }

    string help(){
        return  "Purpose:\n"
                "\tComputes center of mass or geometric center of selection for each frame.\n"
                "\tIf selection is coordinate-dependent updates it every frame.\n"
                "Output:\n"
                "\tFile <label>.dat containing the following columns:\n"
                "\ttime center_x center_y center_z\n"
                "Options:\n"
                "\t--selection <string>\n"
                "\t\tSelection text\n"
                "\t--mass_weighted <true|false>, default: false\n"
                "\t\tCompute center of mass (true) or geometric center (false)";
    }

protected:
    void pre_process(){                
        string sel_text = options("selection").as_string();
        use_mass = options("mass_weighted","false").as_bool();
        sel.modify(system,sel_text);               

        string fname = label+".dat";
        f.open(fname.c_str());
        f << "# time center_x center_y center_z" << endl;
    }

    void process_frame(const Frame_info &info){
        sel.apply();
        f << info.absolute_time << " " << sel.center(use_mass).transpose() << endl;
    }

    void post_process(const Frame_info &info){        
        f.close();
    }    

private:
    Selection sel;
    bool use_mass;
    ofstream f;

};

CREATE_COMPILED_PLUGIN(center)
