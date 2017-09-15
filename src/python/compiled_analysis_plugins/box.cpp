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

TASK_SERIAL(box)
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
        data.clear();
        volume.clear();
    }

    void process_frame(const Frame_info &info) override {
        data.push_back(system.Box(0).extents());
        volume.push_back(system.Box(0).volume());
    }

    void post_process(const Frame_info &info) override {
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
