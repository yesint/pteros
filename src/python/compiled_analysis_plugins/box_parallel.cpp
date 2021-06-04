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
#include "spdlog/fmt/fmt.h"

using namespace std;
using namespace pteros;

struct frame_results: public Result_base {
    Eigen::Vector3f extents;
    float volume;
};

PLUGIN_PARALLEL(box_parallel,frame_results)
public:
    string help(){
        return  "Purpose:\n"
                "\tComputes box vectors and box volume for each frame\n"
                "Output:\n"
                "\tFile <label>.dat containing the following columns:\n"
                "\ttime box_a box_b box_c box_volume\n"
                "Options:\n"
                "\tNone";
    }

protected:
    void pre_process() override {
        //data.clear();
    }

    void process_frame(const Frame_info &info) override {
        //data[info.valid_frame] = frame_results{system.Box(0).extents(),system.Box(0).volume()};
    }

    void post_process(const Frame_info &info) override {        
    }    


};

CREATE_COMPILED_PLUGIN(box_parallel)




