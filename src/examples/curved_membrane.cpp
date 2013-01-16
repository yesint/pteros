#include "pteros/pteros_core.h"
#include "pteros/analysis/trajectory_processor.h"
#include "pteros/analysis/bilayer.h"

#include <unsupported/Eigen/FFT>

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

        for(int i=0;i<po4.size();++i){
            cout << po4.Occupancy(i) << " ";
        }
        cout << endl;
    }

    virtual bool process_frame(const Frame_info& info){
        try{

        cout << "Frame " << info.absolute_frame << endl;
        out.open(string("/media/data/semen/trajectories/grand_challenge/midlines/midline_"
                        + boost::lexical_cast<string>(info.absolute_frame)+".dat").c_str());

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
        float shift = 1.0; // Shift along membrane
        float Ystep = system.Box(0)(1,1)/float(NY+1);
        Vector3f box_dim = system.Box(0).colwise().norm();


        vector<Vector3f> midline; //Membrane midline

        for(;;){
            //cout << iter << endl;
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
            //cout << shift*iter << " " << current_end(0) << " " << current_end(2) << endl;

            if(iter>0){
                Vector3f temp;
                // We need to correct jumps over pbc if they occure
                //cout << "Iter: " << iter << " " << current_end.transpose() << " " << midline.back().transpose() << endl;
                temp = current_end;
                float v;

                v = temp(0)-midline.back()(0);
                if(abs(v) >= box_dim(0)*0.5){
                    (v>0) ? temp(0) -= box_dim(0) : temp(0) += box_dim(0);
                }

                v = temp(2)-midline.back()(1);
                if(abs(v) >= box_dim(2)*0.5){
                    (v>0) ? temp(2) -= box_dim(2) : temp(2) += box_dim(2);
                }

                midline.push_back(Vector3f(temp(0),temp(2),0.0));
            } else {
                midline.push_back(Vector3f(current_end(0),current_end(2),0.0));
            }

            iter++;

            if(iter==1) init_pos = current_end;

            current_end += aver_normal.cross(Y).normalized()*shift;

            if(system.distance(init_pos,current_end,0,true)<shift && iter>10) break;
        }

        // Remove drift of X coordinate
        // We spaned pbc in X direction in iter steps.
        // This corresponds to exactly box_dim(0)
        // We lineraly fit the range (0:box_dim(0)) into iter steps and substract from curve
        float s = midline[0](0);
        for(int i=0;i<midline.size();++i){
            //midline[i](0) += (box_dim(0)-shift)*i/(midline.size()-1);
            midline[i](2) = midline[i](0) - s + (box_dim(0)-shift)*i/(midline.size()-1);
        }

        // Write out
        // Make 20 periodic copies for FFT
        int n = 0;
        for(int c=0;c<1;++c){
            for(int i=0;i<midline.size();++i){
                out << shift*n << " " << midline[i].transpose() << endl;
                ++n;
            }
        }

        cout << "Done in " << iter << " steps" << endl;

        // FFT
        Eigen::FFT<float> fft;
        std::vector<std::complex<float> > freqX, freqZ;
        std::vector<float> dataX, dataZ;
        for(int c=0;c<50;++c){
            for(int i=0;i<midline.size();++i) dataX.push_back(midline[i](2));
            for(int i=0;i<midline.size();++i) dataZ.push_back(midline[i](1));
        }


        fft.fwd(freqX,dataX);
        fft.fwd(freqZ,dataZ);
        ofstream ff(string("/media/data/semen/trajectories/grand_challenge/midlines/FFT_"
                        + boost::lexical_cast<string>(info.absolute_frame)+".dat").c_str());
        for(int i=1;i<freqX.size()/2;++i) ff << 1.0/(shift*i/freqX.size()) << " "
                                             << std::norm(freqX[i]) << " "
                                            << std::norm(freqZ[i]) << endl;
        ff.close();


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

