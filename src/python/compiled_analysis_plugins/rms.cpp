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
                "\t--selection <string>\n"
                "\t\tSelection text"
                "\t--unwrap <float>. Default: 0.2\n"
                "\t\tDo unwrapping of selection based on 'bond distance' criterion"
                "\t\tnegative value means no unwrapping;"
                "\t\tzero means simple nearest neighbour unwrapping,"
                "\t\twhich much faster but fails if selection covers more than 1/2"
                "\t\tof the periodic box size."
                ;
    }

protected:

    void pre_process(){
        mean = 0.0;
        data.clear();
        sel.modify(system, options("selection").as_string() );
        nojump = options("nojump","true").as_bool();
    }

    void process_frame(const pteros::Frame_info &info){
        // Fitting breaks the system, but we have local copy, nobody cares. Cool :)                        

        // Set reference frame for very first processed frame as frame 1        
        if(info.valid_frame==0){
            if(nojump){
                // If nojump is set, do initial unwrapping
                cout << "Initial unwrapping of selection..." << endl;
                float cutoff = 0.2;
                float max_extent = sel.get_system()->Box(0).extents().maxCoeff();
                while(true){
                    try{
                        sel.unwrap_bonds(cutoff);
                    }catch(Pteros_error){
                        cout << "Cutoff " << cutoff << " too small for unwrapping. ";
                        cutoff *= 2.0;
                        cout << "Trying " << cutoff << "..." <<endl;
                        if(cutoff > 0.5*max_extent)
                            throw Pteros_error("Can't unwrap selection with cutoff < 0.5*box_extent.\n"
                                               "Your selection is probably not suitable for RMSD.");
                        continue;
                    }
                    // If we are here unwrapping is successfull
                    break;
                }
                cout << "Unwrapping done." << endl;
            }

            // Create frame 1, which is fixed RMSD reference
            system.frame_dup(0);

            if(nojump){
                // Create frame 2, which is a running reference for unwrapping
                system.frame_dup(0);
            }
        }

        // If nojump is set remove jumps for every atom of selection
        if(nojump){
            for(int i=0;i<sel.size();++i){
                // Get image closest to running reference in frame 2
                sel.XYZ(i,0) = sel.get_system()->
                        Box(0).get_closest_image(sel.XYZ(i,0),sel.XYZ(i,2),false);
                // Update running reference
                sel.XYZ(i,2) = sel.XYZ(i,0);
            }
        }

        // Compute transform with fixed reference in frame 1
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
        f << "# time RMSD:" << endl;
        for(int i=0; i<data.size(); ++i){
            f << i*dt << " " << data[i] << endl;
        }
        f.close();
    }

private:
    std::vector<float> data;
    float mean;
    pteros::Selection sel;
    bool nojump;
};


CREATE_COMPILED_PLUGIN(rms)
