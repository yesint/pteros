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
        // Set beta for all PO4 atoms to 1 - this indicates that
        // All of them are available
        po4.set_beta(0);

        parts.resize(2);                       

        out.open("/media/data/semen/trajectories/grand_challenge/midline.dat");
    }

    virtual bool process_frame(const Frame_info& info){
        try{

        cout << "Frame " << info.absolute_frame << endl;

        Bilayer_point_info inf;
        Vector3f Y(0,1,0);

        if(info.valid_frame ==0 ){
            po4.write("/media/data/semen/trajectories/grand_challenge/sp.pdb");
            cout << "Box:" << system.Box(0) << endl;

            lipids.modify(system,"resname DOPC DOPS");
            bi.create(lipids,"name PO4",2.0);
            inf = bi.point_info(po4.XYZ(0));
            // Determine starting point
            current_end = inf.center;
            cout << "draw sphere \"" << inf.center.transpose()*10 << "\" radius 10" << endl;
            for(int i=0;i<20;++i){
                inf = bi.point_info(current_end);
                current_end = inf.center + inf.normal.cross(Y).normalized()*2.0;
                cout << "draw sphere \"" << current_end.transpose()*10 << "\" radius 10" << endl;
            }

            /*
            cur.modify(system,"name PO4 and beta=0 and distance pbc point "+
                       boost::lexical_cast<string>(current_end(0)) + " " +
                       boost::lexical_cast<string>(current_end(1)) + " " +
                       boost::lexical_cast<string>(current_end(2)) +
                       " < 30"
                      );
            cout << "cur: " << cur.size() << endl;
            cur.set_beta(1); // Exclude
            current_end = cur.center(false,true);
            cout << "Starting end:" << current_end.transpose() << endl;
            out << current_end(0) << " " << current_end(2) << endl;
            */
        }
        // Make selection around current center
        // Include only atoms with beta>0

        /*
        cur.modify(system,"beta=0 and name PO4 and distance pbc point "+
                   boost::lexical_cast<string>(current_end(0)) + " " +
                   boost::lexical_cast<string>(current_end(1)) + " " +
                   boost::lexical_cast<string>(current_end(2)) +
                   " < 31"
                  );

        cout << "cur: " << cur.size() << endl;        
        //cur.set_beta(2);

        // Split obtained selection into two peaces in forward and backward direction
        cur.split_by_connectivity(10,parts);
        //parts[0].set_beta(2);
        //parts[1].set_beta(3);
        //po4.write("/media/data/semen/trajectories/grand_challenge/sp.pdb");

        current_end = parts[0].center(false,true);
        parts[0].set_beta(2);
        inf = bi.point_info(current_end);
        out << inf.center(0) << " " << inf.center(2) << endl;
        cout << "draw sphere \"" << inf.center.transpose()*10 << "\" radius 10" << endl;
        //Vector3f p2 = parts[1].center(false,true);

        //out << p1(0) << " " << p1(2) << endl;
        //out << p2(0) << " " << p2(2) << endl;

        for(int i=0;i<10;++i){
            cur.modify(system,"beta=0 and name PO4 and distance pbc point "+
                       boost::lexical_cast<string>(current_end(0)) + " " +
                       boost::lexical_cast<string>(current_end(1)) + " " +
                       boost::lexical_cast<string>(current_end(2)) +
                       " < 31"
                      );
            cur.set_beta(3+i);
            current_end = cur.center(false,true);
            inf = bi.point_info(current_end);
            out << inf.center(0) << " " << inf.center(2) << endl;
            cout << "draw sphere \"" << inf.center.transpose()*10 << "\" radius 10" << endl;
        }

        po4.write("/media/data/semen/trajectories/grand_challenge/sp.pdb");
        */

        return true;

        } catch(Pteros_error e){
            e.print();
        }
    }

    virtual void post_process(const Frame_info& info){
        out.close();

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

