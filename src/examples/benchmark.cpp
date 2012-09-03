#include "pteros/analysis/trajectory_processor.h"

using namespace std;
using namespace pteros;

// RMSD fitting
class Bench1: public Consumer {
public:
    Bench1(Trajectory_processor* pr): Consumer(pr){
    }
protected:
    virtual void pre_process(){
        sel.modify(system, "all");
    }

    virtual bool process_frame(const Frame_info& info){
        if(info.valid_frame==0){
            system.frame_dup(0);
        }

        Eigen::Affine3f t = sel.fit_transform(0,1);
        sel.apply_transform(t);
        cout << "RMSD of frame " << info.valid_frame << " " << sel.rmsd(0,1) << endl;
        return true;
    }

    Selection sel;
};

// Contacts computation
class Bench2: public Consumer {
public:
    Bench2(Trajectory_processor* pr): Consumer(pr){
    }
protected:
    virtual void pre_process(){
        sel1.modify(system, "resid 1 to 100");
        sel2.modify(system, "resid 102 to 200");
    }

    virtual bool process_frame(const Frame_info& info){
        vector<Eigen::Vector2i> bon;
        bon.clear();
        Grid_searcher(0.5,sel1,sel2,bon,true);

        cout << "frame " << info.valid_frame << " bonds: " << bon.size() << endl;
        return true;
    }

    Selection sel1, sel2;
};


// Selecting each residue
class Bench3: public Consumer {
public:
    Bench3(Trajectory_processor* pr): Consumer(pr){
    }
protected:
    virtual void pre_process(){

    }

    virtual bool process_frame(const Frame_info& info){
        vector<Selection> sel;
        Selection all(system,"all");
        all.each_residue(sel);

        cout << "frame " << info.valid_frame << endl;
        return true;
    }


};


int main(int argc, char** argv){

    try{


    Options_tree options;
    options.from_command_line(argc,argv);
    int num = options.get_value<int>("bench");
    cout << num << endl;
    if(num>3){
        System s("/home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_pdb2gmx.gro");
        Selection sel(s,"all");

        std::clock_t start;
        double duration;
        start = std::clock();


        Eigen::Vector3f v(0.1,0.1,0.1);
        Eigen::Vector3f m1,m2;

        for(int i=0; i<100000; ++i){
            switch(num){
            case 4: {
                sel.translate(v);
                break;
            }
            case 5: {
                sel.rotate(0,0.1);
                break;
            }
            case 6: {
                sel.center();
                break;
            }
            case 7: {
                sel.minmax(m2,m2);
                break;
            }
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;

        return 1;
    }

    Trajectory_processor engine(options);

    boost::shared_ptr<Consumer> b1,b2,b3;

    switch(num){
    case 1: {
        b1 = boost::shared_ptr<Bench1>(new Bench1(&engine));
        break;
    }
    case 2: {
        b2 = boost::shared_ptr<Bench2>(new Bench2(&engine));
        break;
    }
    case 3: {
        b3 = boost::shared_ptr<Bench3>(new Bench3(&engine));
        break;
    }
    case 0: {
        b1 = boost::shared_ptr<Bench1>(new Bench1(&engine));
        b2 = boost::shared_ptr<Bench2>(new Bench2(&engine));
        b3 = boost::shared_ptr<Bench3>(new Bench3(&engine));
        break;
    }
    }

    engine.run();

    } catch(Pteros_error e){
        cout << e.what() << endl;
    }

}

