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
                "\t\tSelection text\n"
                "\t-nojump <true|false>. Default: true\n"
                "\t\tRemove jumps of atoms over periodic box boundary.\n"
                "\t\tSetting this to false is only meaningful is very special cases."
                ;
    }

protected:

    void pre_process(){        
        data.clear();        

        // rms_sel is required
        rms_sel.modify(system, options("rms_sel").as_string() );

        // fit_sel is optional
        string fit_str = options("fit_sel","").as_string();
        if(fit_str!="") fit_sel.modify(system, fit_str);

        if((fit_sel.size()>0 && fit_sel.size()<3) || (fit_sel.size()==0 && rms_sel.size()<3)){
            throw Pteros_error("Can't fit selection with less than 3 atoms!");
        }

        // Add our selections to nojump list
        add_no_jump_atoms(fit_sel);
        add_no_jump_atoms(rms_sel);

        // Create frame 1 for fitting
        system.frame_dup(0);
    }     

    void process_frame(const pteros::Frame_info &info){                

        // Compute RMSD with fixed reference in frame 1
        if(fit_sel.size()>2){
            Eigen::Affine3f trans = fit_sel.fit_transform(0,1);
            rms_sel.apply_transform(trans);
        } else {
            rms_sel.fit(0,1);
        }

        float v = rms_sel.rmsd(0,1);
        data.push_back(v);        
    }

    void post_process(const pteros::Frame_info &info){        
        // Output
        string fname = label+".dat";
        // Get time step in frames and time
        float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);

        ofstream f(fname.c_str());
        f << "# RMSD of selection [" << rms_sel.get_text() << "]" << endl;
        if(fit_sel.size()>2){
            f << "# after fitting of selection [" << fit_sel.get_text() << "]" << endl;
        }
        f << "# time(ns) RMSD(nm)" << endl;
        for(int i=0; i<data.size(); ++i){
            f << i*dt << " " << data[i] << endl;
        }
        f.close();
    }

private:
    vector<float> data;
    Selection fit_sel, rms_sel;
};


CREATE_COMPILED_PLUGIN(rms)
