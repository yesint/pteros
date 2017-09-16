/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
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
