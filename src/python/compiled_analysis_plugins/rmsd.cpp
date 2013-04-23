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

class rmsd: public pteros::Compiled_plugin_base {
public:
    rmsd(pteros::Trajectory_processor* pr, pteros::Options_tree* opt): Compiled_plugin_base(pr,opt) {}
protected:

    void pre_process(){
        mean = 0.0;
        data.clear();
        sel.modify(system, options->get_value<string>("selection") );
    }

    void process_frame(const pteros::Frame_info &info){
        // Fitting breaks the system, but we have local copy, nobody cares. Cool :)
        // Set reference frame
        if(info.valid_frame==0){
            system.frame_dup(0);
        }

        Eigen::Affine3f trans = sel.fit_transform(0,1);
        sel.apply_transform(trans);
        float v = sel.rmsd(0,1);
        data.push_back(v);
        mean += v;
    }

    void post_process(const pteros::Frame_info &info){
        mean /= (float)info.valid_frame;
        // Output
        string fname = label+".dat";
        // Get time step in frames and time
        float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);

        ofstream f(fname.c_str());
        f << "# RMSD of selection [" << sel.get_text() << "]" << endl;
        f << "# Mean: " << mean << endl;
        f << "# time X Y Z:" << endl;
        for(int i=0; i<data.size(); ++i){
            f << i*dt << " " << data[i] << endl;
        }
        f.close();
    }

private:
    std::vector<float> data;
    float mean;
    pteros::Selection sel;
};


CREATE_COMPILED_PLUGIN(rmsd)
