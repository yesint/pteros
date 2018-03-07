#include "pteros/core/pteros_error.h"
#include "pteros/analysis/options.h"
#include "pteros/pteros.h"
#include <Eigen/Core>
#include <fstream>

using namespace std;
using namespace pteros;
using namespace Eigen;

int main(int argc, char** argv){
    try {
    Options opt;
    parse_command_line(argc,argv,opt);

    string fname;

    // Load structure
    fname = opt("struct").as_string();
    System sys(fname);
    // Get indexes of ROH atoms of CHOL molecules
    Selection ROH(sys,std::string("name ROH"));
    vector<int> ROH_index = ROH.get_index();

    // Occupancy of monolayers over time averaged over all trajectories
    vector<int> mon1_occ, mon2_occ;

    // Cycle over supplied trajectories
    vector<string> traj_list = opt("traj").as_strings();
    for(string& ff: traj_list){
        // Delete all frames
        sys.frame_delete();
        // Load trajectory
        sys.load(ff);
        // Occupancy of monolayers for this trajectory over time
        vector<int> occ1, occ2;
        // Cycle over frames
        Selection sel(sys);
        string sel_text;
        for(int fr=0; fr<sys.num_frames(); ++fr){
            cout << "Frame " << fr << endl;
            occ1.push_back(0);
            occ2.push_back(0);
            // For each lipid find local center of bilayer
            for(int lip=0; lip<ROH.size(); ++lip){
                sel_text = "not resname W and within_xy 1.5 of index " + to_string(ROH_index[lip]);
                sel.modify(sel_text);
                sel.set_frame(fr);
                float localZ = sel.center()(2);
                // See if current ROH is above or below localZ
                if(ROH.z(lip,fr)<localZ){
                    ++occ1.back();
                } else {
                    ++occ2.back();
                }
            } //Over lipids
        } // Over frames

        if(mon1_occ.size()==0) mon1_occ.resize(occ1.size());
        if(mon2_occ.size()==0) mon2_occ.resize(occ2.size());
        // Add to global occupancies
        for(int i=0;i<mon1_occ.size();++i){
            mon1_occ[i] += occ1[i];
            mon2_occ[i] += occ2[i];
        }

    } // Over trajectories

    // Average occupancies
    float N = traj_list.size();

    // Write occupancies
    ofstream f("in_time.dat");
    for(int i=0;i<mon1_occ.size();++i){
        f << i << " " << float(mon1_occ[i])/N << " " << float(mon2_occ[i])/N << " "
          << float(mon1_occ[i])/float(mon2_occ[i]) << endl;
    }
    f.close();

    } catch(const Pteros_error& e){
        cout << e.what() << endl;
    }
}

