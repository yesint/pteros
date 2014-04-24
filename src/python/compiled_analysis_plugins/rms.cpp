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
#include "pteros/core/pteros_error.h"
#include <fstream>

using namespace std;
using namespace pteros;

class rms: public pteros::Compiled_plugin_base {
public:
    rms(pteros::Trajectory_processor* pr, const pteros::Options& opt): Compiled_plugin_base(pr,opt) {}

    string help(){
        return  "Purpose:\n"
                "\tComputes RMSD of each frame for given selection.\n"
                "\tThe first loaded frame is used as a reference.\n"
                "\tSelection should be coordinate-independent.\n"
                "Output:\n"
                "\tFile <label>.dat containing the following columns:\n"
                "\ttime RMSD\n"
                "\tAlso reports mean RMSD in the file header.\n"
                "Options:\n"
                "\t-selection <string>\n"
                "\t\tSelection text"
                "\t-nojump <true|false>. Default: true\n"
                "\t\tRemove jumps of atoms aver periodic box boundary."
                "\t\tSetting this to false is only meaningful is very special cases."
                ;
    }

protected:

    void pre_process(){        
        data.clear();
        sel.modify(system, options("selection").as_string() );
        nojump = options("nojump","true").as_bool();

        // Add our selection to nojump list if asked
        if(nojump){
            add_no_jump_atoms(sel);
        }
    }     

    void process_frame(const pteros::Frame_info &info){                
        // Fitting breaks the system, but we have local copy, nobody cares. Cool :)                        
        if(info.valid_frame==0){
            // Create frame 1 for fitting
            system.frame_dup(0);
        }

        // Compute transform with fixed reference in frame 1
        Eigen::Affine3f trans = sel.fit_transform(0,1);
        sel.apply_transform(trans);
        float v = sel.rmsd(0,1);

        data.push_back(v);        
    }

    void post_process(const pteros::Frame_info &info){        
        // Output
        string fname = label+".dat";
        // Get time step in frames and time
        float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);

        ofstream f(fname.c_str());
        f << "# RMSD of selection [" << sel.get_text() << "]" << endl;        
        f << "# time RMSD:" << endl;
        for(int i=0; i<data.size(); ++i){
            f << i*dt << " " << data[i] << endl;
        }
        f.close();
    }

private:
    std::vector<float> data;    
    pteros::Selection sel;
    bool nojump;
};


CREATE_COMPILED_PLUGIN(rms)
