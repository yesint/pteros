/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#include "pteros/python/compiled_plugin.h"
#include <fstream>
#include "fmt/core.h"

using namespace std;
using namespace pteros;

TASK_SERIAL(box)
public:
    string help() override {
        return
R"(Purpose:
    Computes box vectors and box volume for each frame
Output:
    File <label>.dat containing the following columns:
    time box_a box_b box_c box_volume
Options:
    None
)";
    }

protected:
    void pre_process() override {
        data.clear();
        volume.clear();
    }

    void process_frame(const FrameInfo &info) override {
        data.push_back(system.box(0).extents());
        volume.push_back(system.box(0).volume());
    }

    void post_process(const FrameInfo &info) override {
        // Output
        string fname = fmt::format("box_id{}.dat",get_id());
        // Get time step in frames and time
        float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);

        ofstream f(fname.c_str());
        f << "# time box_a box_b box_c box_volume:" << endl;
        for(int i=0; i<data.size(); ++i){
            f << i*dt << " " << data[i].transpose() << " " << volume[i] << endl;
        }
        f.close();
    }    

private:
    vector<Eigen::Vector3f> data;
    vector<float> volume;
};

CREATE_COMPILED_PLUGIN(box)




