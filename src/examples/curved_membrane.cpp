#include "pteros/pteros_core.h"
#include "pteros/analysis/trajectory_processor.h"
#include "pteros/analysis/bilayer.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


class Membrane_processor: public Consumer {
public:
    Membrane_processor(Trajectory_processor* pr, Options_tree* opt): Consumer(pr) {
        //options = opt;
    }
protected:
    virtual void pre_process(){
        po4.modify(system,"name PO4");

        cout << "PO4: " << po4.size() << endl;        
    }

    virtual bool process_frame(const Frame_info& info){
        try{

        cout << "Frame " << info.valid_frame << endl;
        out.open(string("/media/data/semen/trajectories/grand_challenge/midline_"
                        + boost::lexical_cast<string>(info.valid_frame)+".dat").c_str());

        Bilayer_point_info inf;
        Vector3f Y(0,1,0);
        Vector3f init_pos;
        Vector3f tmp;

        if(info.valid_frame ==0 ){
            po4.write("/media/data/semen/trajectories/grand_challenge/sp.pdb");
            cout << "Box:" << system.Box(0) << endl;

            lipids.modify(system,"resname DOPC DOPS");
            bi.create(lipids,"name PO4",2.0);
        }

        inf = bi.point_info(po4.XYZ(0));
        // Determine starting point
        current_end = inf.center;
        //cout << "draw sphere \"" << inf.center.transpose()*10 << "\" radius 10" << endl;

        int iter = 0;
        int NY = 5;
        float shift = 2.0;
        float Ystep = system.Box(0)(1,1)/float(NY+1);
        for(;;){
            iter++;
            // We are computing average of bilayer centers by going along Y
            // Taking 10 points distributed over whole Y dimension
            tmp = current_end;
            tmp(1) = 0;
            current_end.fill(0);
            Vector3f aver_normal;
            aver_normal.fill(0);
            for(int i=0;i<NY;++i){
                inf = bi.point_info(tmp);
                current_end += inf.center;
                aver_normal += inf.normal;
                tmp(1) += Ystep;
            }
            current_end /= float(NY);
            //cout << "draw sphere \"" << current_end.transpose()*10 << "\" radius 10" << endl;
            out << current_end(0) << " " << current_end(2) << endl;

            if(iter==1) init_pos = current_end;

            current_end += aver_normal.cross(Y).normalized()*shift;

            if((init_pos-current_end).norm()<shift && iter>10) break;
        }
        cout << "Done in " << iter << " steps" << endl;
        out.close();

        return true;

        } catch(Pteros_error e){
            e.print();
        }
    }

    virtual void post_process(const Frame_info& info){

    }

    Vector3f current_end;
    Selection po4;
    Selection cur;
    Selection lipids;
    float dl;
    vector<float> midline_x;
    vector<float> midline_z;
    vector<Selection> parts;
    Bilayer bi;
    ofstream out;
};


int main(int argc, char** argv){

    try {
        Options_tree options;
        options.from_command_line(argc,argv);
        Trajectory_processor proc(options);
        Membrane_processor membr(&proc,&options);
        proc.run();
    } catch(Pteros_error e){
        e.print();
    }

}

