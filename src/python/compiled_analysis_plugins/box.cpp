/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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

using namespace std;
using namespace pteros;

class box: public Compiled_plugin_base {
public:
    box(Trajectory_processor* pr, Options_tree* opt): Compiled_plugin_base(pr,opt) {}

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
    void pre_process(){
        data.clear();
        volume.clear();
    }

    void process_frame(const Frame_info &info){
        Eigen::Matrix3f m = system.Box(0);
        data.push_back(m.diagonal());
        volume.push_back(m.diagonal().prod());        
    }

    void post_process(const Frame_info &info){
        // Output
        string fname = label+".dat";
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
