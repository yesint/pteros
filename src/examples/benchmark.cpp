#include "pteros/core/distance_search.h"
#include "pteros/core/pteros_error.h"
#include "pteros/analysis/task_plugin.h"
#include "pteros/analysis/trajectory_reader.h"
#include <ctime>

using namespace std;
using namespace pteros;

// RMSD fitting

TASK_SERIAL(Bench1)
protected:
    virtual void pre_process() override {
        sel.modify(system,std::string("all"));
    }

    virtual void process_frame(const Frame_info& info) override {
        if(info.valid_frame==0){
            system.frame_dup(0);
        }

        Eigen::Affine3f t = sel.fit_transform(0,1);
        sel.apply_transform(t);
        cout << "RMSD of frame " << info.valid_frame << " " << sel.rmsd(0,1) << endl;        
    }

    virtual void post_process(const Frame_info& info) override {}

    Selection sel;
};

// Contacts computation
TASK_SERIAL(Bench2)
protected:
    virtual void pre_process() override {
        sel1.modify(system,std::string("resid 1 to 100"));
        sel2.modify(system,std::string("resid 102 to 200"));
    }

    virtual void process_frame(const Frame_info& info) override {
        vector<Eigen::Vector2i> bon;
        bon.clear();
        search_contacts(0.5,sel1,sel2,bon,true);

        cout << "frame " << info.valid_frame << " bonds: " << bon.size() << endl;        
    }

    virtual void post_process(const Frame_info& info) override {}


    Selection sel1, sel2;
};


// Selecting each residue
TASK_SERIAL(Bench3)
protected:
    virtual void pre_process() override  {

    }

    virtual void process_frame(const Frame_info& info) override {
        vector<Selection> sel;
        Selection all(system,std::string("all"));
        all.each_residue(sel);

        cout << "frame " << info.valid_frame << endl;        
    }

    virtual void post_process(const Frame_info& info) override {}
};


int main(int argc, char** argv){

    try{

    Options options;
    parse_command_line(argc,argv,options);
    int num = options("bench").as_int();
    cout << num << endl;
    if(num>3){
        System s("/home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_pdb2gmx.gro");
        Selection sel(s,std::string("all"));

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
            case 5:
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

    Trajectory_reader engine(options);





    switch(num){
    case 1: {
        engine.add_task( new Bench1(options) );
        break;
    }
    case 2: {
        engine.add_task( new Bench2(options) );
        break;
    }
    case 3: {
        engine.add_task( new Bench3(options) );
        break;
    }
    case 0: {
        engine.add_task( new Bench1(options) );
        engine.add_task( new Bench2(options) );
        engine.add_task( new Bench3(options) );
        break;
    }
    }

    engine.run();

    } catch(const Pteros_error& e){
        cout << e.what() << endl;
    }

}

