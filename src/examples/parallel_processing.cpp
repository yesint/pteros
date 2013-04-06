#include "pteros/analysis/trajectory_processor.h"
#include "pteros/analysis/consumer.h"
#include "pteros/core/grid_search.h"

using namespace std;
using namespace pteros;

class Our_task: public Consumer {
public:
    // Constructor
    Our_task(Trajectory_processor* pr, string s1, string s2, float d): Consumer(pr){
        sel1_text = s1;
        sel2_text = s2;
        dist = d;
    }
protected:      
    virtual void pre_process(){
        sel1.modify(system,sel1_text);
        sel2.modify(system,sel2_text);
    }

    virtual bool process_frame(const Frame_info& info){
        vector<Eigen::Vector2i> contacts;
        Grid_searcher(dist,sel1,sel2,contacts);
        float en = system.non_bond_energy(contacts,0).total;
        cout << "Energy for time" << info.absolute_time << " is " << en << endl;
        return true;
    }

    // Variables, which are specific for our analysis
    Selection sel1, sel2;
    string sel1_text, sel2_text;
    float dist;
};

int main(int argc, char** argv){
    Options_tree options;
    options.from_command_line(argc,argv);
    Trajectory_processor engine(options);
    string s1 = options.get_value<string>("selection1");
    string s2 = options.get_value<string>("selection2");
    float d = options.get_value<float>("distance");
    Our_task task(&engine,s1,s2,d);
    engine.run();
}

