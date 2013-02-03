#include "pteros/analysis/trajectory_processor.h"
#include "pteros/analysis/bilayer.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

// Data for each time frame
struct Chol_data {
    float t;
    int n1,n2;
    float flip_time;
    float mean_curvature; // abs computed for both parts
    vector<int> domain; // domain for each chol
    vector<int> monolayer; // monolayer for each chol
    vector<float> center_dist; // Distance from center for each chol
    vector<float> curvature; // Curvature for each chol

    string print(){
        return boost::lexical_cast<string>(t) + " "
                + boost::lexical_cast<string>(n1) + " "
                + boost::lexical_cast<string>(n2) + " "
                + boost::lexical_cast<string>(flip_time) + " "                
                + boost::lexical_cast<string>(mean_curvature) + " "
                ;
    }
};


class Chol_counter: public Consumer {
public:
    Chol_counter(Trajectory_processor* pr, Options_tree* opt): Consumer(pr) {
        options = opt;
    }
protected:
    virtual void pre_process(){
        bilayer.modify(system,options->get_value<string>("lipids_selection"));
        roh.modify(system,options->get_value<string>("chol_head_selection"));
        lip_name1 = options->get_value<string>("lipid_name1");
        lip_name2 = options->get_value<string>("lipid_name2");

        file_path = options->get_value<string>("file_path",".");

        // Parse bilayer
        marker_sel_text = options->get_value<string>("lipid_marker_selection");
        bi.create(bilayer,marker_sel_text,2.0);        

        out.open(options->get_value<string>("output_file","trace.dat").c_str());
    }

    virtual bool process_frame(const Frame_info& info){
        // Curvature values for each chol and separately for each lipid domain type        

        vector<float> curv_hist_data, curv_hist_data1, curv_hist_data2;


        bilayer.set_frame(0);

        Chol_data dum;
        trace.push_back(dum);
        trace.back().domain.resize(roh.size());
        //trace.back().num_type1.resize(roh.size());
        //trace.back().num_type2.resize(roh.size());
        trace.back().monolayer.resize(roh.size());
        trace.back().center_dist.resize(roh.size());
        trace.back().curvature.resize(roh.size());
        trace.back().t = info.absolute_time;

        // Read midline
        vector<float> curv;
        vector<Vector2f> midline;

        float tmp;

        // Read midline from file
        ifstream ff(string(file_path
                           +"/midline_"
                           +boost::lexical_cast<string>(info.absolute_frame)
                           +".dat").c_str());

        if(!ff.good()){
            cout << "Error reading file!" << endl;
            ff.close();
            return true;
        }

        string line;
        stringstream ss;
        float x,z,c;
        while(getline(ff,line)){
            ss.clear();
            ss.str(line);
            ss >> tmp >> x >> z >> tmp >> c;
            midline.push_back(Vector2f(x,z));
            curv.push_back(c);
        }
        ff.close();


        // Cycle over all roh atoms and compute their positions
        int n1=0, n2=0, Nflip=0, non_flip=0;
        for(int i=0; i<roh.size(); ++i){

            // Get info about current position
            Bilayer_point_info data = bi.point_info(roh.XYZ(i));            

            // Analyze which lipids surround this chol
            boost::shared_ptr<Selection> ptr;
            int num1=0, num2=0;

            if(data.monolayer==1){                
                ptr = data.spot1_ptr;                
            } else {                
                ptr = data.spot2_ptr;
            }

            for(int j=0; j<ptr->size(); ++j){
                if(ptr->Resname(j)==lip_name1) ++num1;
                if(ptr->Resname(j)==lip_name2) ++num2;
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

            // Curvature accounting

            // Get closest midline point for current bilayer center
            float d,min_d = 1e20;
            int min_ind = -1;
            Vector3f point;
            for(int j=0;j<midline.size();++j){
                point(0) = midline[j](0);
                point(1) = data.center(1);
                point(2) = midline[j](1);

                d = system.distance(point,data.center,0,true);

                if(d<min_d){
                    min_d = d;
                    min_ind = j;
                }
            }
            // Store curvature
            trace.back().curvature[i] = curv[min_ind];

            //cout << "chol #" << i << " " << min_ind << " " << curv[min_ind] << " " << curv.size() << endl;

            // Store histograms data
            /*
            if(min_ind<midline.size()/2){
                tmp = curv[min_ind];
            } else {
                tmp = -curv[min_ind];
            }
            */
            tmp = abs(curv[min_ind]);

            curv_hist_data.push_back(tmp);
            if(trace.back().domain[i] == 1) curv_hist_data1.push_back(tmp);
            if(trace.back().domain[i] == 2) curv_hist_data2.push_back(tmp);
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

        // Compute mean curvature. Two parts must be accounted with opposite signs
        tmp = 0.0;
        for(int i=0;i<midline.size();++i){
            //(i<midline.size()/2) ? tmp += curv[i] : tmp -= curv[i];
            tmp += curv[i]*curv[i];
        }
        tmp = sqrt(tmp/float(midline.size()));
        trace.back().mean_curvature = tmp;

        cout << "Frame " << info.absolute_time << " " << trace.back().print() << endl;

        /*
        // Output
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
        */
        out << trace.back().print() << endl;
        //ff << " " << mean1/m1 << " " << mean2/m2 << endl;

        // Output histograms
        ofstream h(string(file_path
                          +"/hist_all_"
                          +boost::lexical_cast<string>(info.absolute_frame)
                          +".dat").c_str());
        for(int j=-0;j<curv_hist_data.size();++j) h << curv_hist_data[j] << endl;
        h.close();

        h.open(string(file_path
                          +"/hist_1_"
                          +boost::lexical_cast<string>(info.absolute_frame)
                          +".dat").c_str());
        for(int j=-0;j<curv_hist_data1.size();++j) h << curv_hist_data1[j] << endl;
        h.close();

        h.open(string(file_path
                          +"/hist_2_"
                          +boost::lexical_cast<string>(info.absolute_frame)
                          +".dat").c_str());
        for(int j=-0;j<curv_hist_data2.size();++j) h << curv_hist_data2[j] << endl;
        h.close();

        return true;
    }

    virtual void post_process(const Frame_info& info){

        // Write trace        
        /*
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
        */
        out.close();
    }

    Bilayer bi;
    Selection bilayer; //Bilayer
    Selection roh; // Cholesterol heads
    string lip_name1, lip_name2; // Name of the lipids

    string marker_sel_text;
    Options_tree* options;
    // Trace in time
    vector<Chol_data> trace;    
    float last_t;

    string file_path;
    ofstream out;
};

int main(int argc, char** argv){
    try {
        Options_tree options;
        options.from_command_line(argc,argv);
        Trajectory_processor proc(options);
        Chol_counter counter(&proc,&options);
        proc.run();
    } catch(Pteros_error e){
        e.print();
    }
}

