#include "pteros/analysis/compiled_analysis_task.h"
#include <fstream>

using namespace std;
using namespace pteros;

PLUGIN_SERIAL(center)
public:

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
        use_mass = options("mass","false").as_bool();
        string fname = get_id()+".dat";
        f.open(fname.c_str());
        f << "# time center_x center_y center_z" << endl;
        string sel_text = options("sel").as_string();
        sel.modify(system,sel_text);        
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

COMPILED_ANALYSIS_TASK(center,nullptr)

