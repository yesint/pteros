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

//COMPILED_ANALYSIS_TASK(center,nullptr)

int main(int argc, char** argv){
    try {
        Options options;
        parse_command_line(argc,argv,options);
        Trajectory_reader engine(options);
        auto task = new center(options);
        engine.add_task(task);
        //engine.register_collector(_collector);
        cout << "-------------------------------------------------------------" << endl;
        cout << "  This is stand-alone Pteros analysis plugin " << endl;
        cout << "-------------------------------------------------------------" << endl;
        if(!options.has("f") && !options.has("help")){
            cout << "Usage:" << endl;
            cout << "\n\tFor specific task options use '-help task'" << endl;
            cout << "\tFor trajectory processing options use '-help traj'" << endl;
            cout << "\tFor all available options use '-help all' or just '-help'" << endl;
            return 1;
        }
        if(options.has("help")){
            string help = options("help","").as_string();
            if(help=="traj"){
                cout << engine.help() << endl;
            } else if(help=="task"){
                cout << task->help() << endl;
            } else {
                cout << task->help() << endl << endl;
                cout << engine.help() << endl;
            }
            return 1;
        }

        engine.run();
    } catch(const Pteros_error& e) {
        LOG()->error(e.what());
    } catch (const std::exception& e) {
        LOG()->error(e.what());
    } catch(...) {
        LOG()->error("Unknown error");
    }
}
