#include "bilayer.h"
#include "pteros/core/pteros_error.h"
#include "pteros/pteros.h"
#include "pteros/analysis/trajectory_reader.h"
#include "pteros/analysis/task_plugin.h"
#include <fstream>

using namespace std;
using namespace pteros;
using namespace Eigen;

// Data for each time frame
struct Chol_data {
    float t;
    int n1,n2;
    float flip_time;
    vector<int> domain; // domain for each chol
    vector<int> monolayer; // monolayer for each chol
    vector<float> center_dist; // Distance from center for each chol    

    string print(){
        return    to_string(t) + " "
                + to_string(n1) + " "
                + to_string(n2) + " "
                + to_string(flip_time) + " "
                ;
    }
};


TASK_SERIAL(Chol_counter)
protected:
    virtual void pre_process() override {
        bilayer.modify(system,options("lipids_selection").as_string());
        roh.modify(system,options("chol_head_selection").as_string());
        lip_name1 = options("lipid_name1").as_string();
        lip_name2 = options("lipid_name2").as_string();

        // Parse bilayer
        marker_sel_text = options("lipid_marker_selection").as_string();
        bi.create(bilayer,marker_sel_text,2.0);        
    }

    virtual void process_frame(const Frame_info& info) override {
        bilayer.set_frame(0);

        Chol_data dum;
        trace.push_back(dum);
        trace.back().domain.resize(roh.size());
        //trace.back().num_type1.resize(roh.size());
        //trace.back().num_type2.resize(roh.size());
        trace.back().monolayer.resize(roh.size());
        trace.back().center_dist.resize(roh.size());
        trace.back().t = info.absolute_time;

        // Cycle over all roh atoms and compute their positions
        int n1=0, n2=0, Nflip=0, non_flip=0;
        for(int i=0; i<roh.size(); ++i){
            // Get info about current position
            Bilayer_point_info data = bi.point_info(roh.xyz(i));

            // Analyze which lipids surround this chol
            std::shared_ptr<Selection> ptr;
            int num1=0, num2=0;

            if(data.monolayer==1){                
                ptr = data.spot1_ptr;                
            } else {                
                ptr = data.spot2_ptr;
            }

            for(int j=0; j<ptr->size(); ++j){
                if(ptr->resname(j)==lip_name1) ++num1;
                if(ptr->resname(j)==lip_name2) ++num2;
            }

            if(num1>num2){                
                trace.back().domain[i] = 1;
                ++n1;
            } else if(num1<num2){
                trace.back().domain[i] = 2;
                ++n2;
            } else {
                trace.back().domain[i] = 0; // Undefined = on the domain boundary
            }

            // Track monolayer
            trace.back().monolayer[i] = data.monolayer;
            // Track center dist
            trace.back().center_dist[i] = data.center_dist;
            
            // Flip-flop accounting
            if(info.valid_frame>0){
                // Can't do this on first frame
                int prev = trace.size()-2; // Previoud frame
                int cur = trace.size()-1; // This frame

                bool in_same_domain = (trace[prev].domain[i] == trace[cur].domain[i]);
                bool in_same_monolayer = (trace[prev].monolayer[i] == trace[cur].monolayer[i]);
                bool any_on_boundary = (trace[prev].domain[i]*trace[cur].domain[i] == 0);

                if(!in_same_domain && !in_same_monolayer && !any_on_boundary){
                    ++Nflip;
                } else if(in_same_domain && in_same_monolayer && !any_on_boundary) {
                    ++non_flip;
                }
            }

        }

        if(info.valid_frame>0){
            // Compute estimated flip-flop time
            float dt = info.absolute_time-last_t;
            //cout << possible_flip << " " << Nflip << " " << dt << endl;

            //trace.back().flip_time = (float)(non_flip+Nflip)*dt/(float)Nflip;
            trace.back().flip_time = (float)Nflip/(float)(non_flip+Nflip);
            //cout << (float)(non_flip+Nflip)/(float)Nflip << " --- " << info.absolute_time << " " << last_t << endl;
        } else {
            trace.back().flip_time = -1; //Invalid on first frame
        }

        last_t = info.absolute_time;

        trace.back().n1 = n1;
        trace.back().n2 = n2;

        cout << "Frame " << info.absolute_frame << " " << trace.back().print() << endl;
    }

    virtual void post_process(const Frame_info& info) override {

        // Write trace
        ofstream ff(options("output_file","trace.dat").as_string().c_str());
        for(int i=0; i<trace.size();++i){
            ff << trace[i].print();
            // Compute mean center distance for each domain
            float mean1 = 0.0, mean2 = 0.0;
            float m1=0,m2=0;
            for(int j=0; j<roh.size(); ++j){
                if(trace[i].domain[j]==1){
                    mean1+=trace[i].center_dist[j];
                    ++m1;
                }
                if(trace[i].domain[j]==2){
                    mean2+=trace[i].center_dist[j];
                    ++m2;
                }
            }
            ff << " " << mean1/m1 << " " << mean2/m2 << endl;
        }               
        ff.close();        
    }

    Bilayer bi;
    Selection bilayer; //Bilayer
    Selection roh; // Cholesterol heads
    string lip_name1, lip_name2; // Name of the lipids

    string marker_sel_text;
    Options options;
    // Trace in time
    vector<Chol_data> trace;    
    float last_t;
};

int main(int argc, char** argv){
    try {
        Options options;
        parse_command_line(argc,argv,options);
        Trajectory_reader proc(options);
        proc.add_task(new Chol_counter(options));
        proc.run();
    } catch(const Pteros_error& e){
        cout << e.what() << endl;
    }
}

