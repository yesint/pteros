/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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

TASK_SERIAL(rms)
public:    

    string help() override {
        return  "Purpose:\n"
                "\tComputes RMSD of each frame for given selection.\n"                
                "\tSelection should be coordinate-independent.\n"
                "Output:\n"
                "\tFile rms_id<id>.dat containing the following columns:\n"
                "\ttime RMSD\n"
                "\tAlso reports mean RMSD in the file header.\n"
                "Options:\n"
                "\t-fit_sel <string>\n"
                "\t\tFitting selection text\n"

                "\t-rms_sel <string string...>\n"
                "\t\tRMSD selections texts\n"

                "\t-nojump <distance>. Default: 0\n"
                "\t\tRemove jumps of atoms over periodic box boundary.\n"
                "\t\tAtoms, which should not jump, are unwrapped with\n"
                "\t\tgiven distance on the first frame.\n"
                "\t\tZero means find unwrap distance automatically.\n"
                "\t\tDistance -1 means no jump removal\n"
                "\t\tthis is only meaningful is very special cases!\n"
                ;
    }

protected:

    void pre_process() override {
        data.clear();        

        // rms_selections
        vector<string> strs = options("rms_sel").as_strings();
        if(strs.size()==0) throw Pteros_error("At least one rms selection required!");
        // Create selections
        for(auto& s: strs) rms_sel.emplace_back(system,s);

        if(check_selection_overlap(rms_sel)) throw Pteros_error("Selections should not overlap!");

        // fit_sel is optional
        string fit_str = options("fit_sel","").as_string();
        if(fit_str!=""){
            fit_sel.modify(system, fit_str);
        } else {
            log->info("Using first rms selection for fitting");
            fit_sel = rms_sel[0];
        }

        if(fit_sel.size()<3){
            throw Pteros_error("Can't fit selection with less than 3 atoms!");
        }

        float d = options("nojump","0").as_float();
        if(d>=0){
            // Add our selections to nojump list
            jump_remover.add_atoms(fit_sel);
            for(auto& sel: rms_sel) jump_remover.add_atoms(sel);
            jump_remover.set_unwrap_dist(d);
        }

        data.resize(rms_sel.size());
    }     

    void process_frame(const pteros::Frame_info &info) override {
        if(info.valid_frame==0){
            // Create frame 1 for fitting
            system.frame_dup(0);
        }

        // Compute RMSD with fixed reference in frame 1        
        auto trans = fit_sel.fit_transform(0,1);
        for(int i=0; i<rms_sel.size(); ++i){
            rms_sel[i].apply_transform(trans);
            float v = rms_sel[i].rmsd(0,1);
            data[i].push_back(v);
        }
    }

    void post_process(const pteros::Frame_info &info) override {
        // Output
        string fname = fmt::format("rms_id{}.dat",get_id());
        // Get time step in frames and time
        float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);

        ofstream f(fname.c_str());
        f << "# RMSD of selections:"<<endl;
        for(int i=0; i<rms_sel.size(); ++i){
            f << "# " << i << ": '" << rms_sel[i].get_text().substr(0,80) << "'" << endl;
        }
        f << "# after fitting of selection '" << fit_sel.get_text() << "'" << endl;
        f << "# time(ns) RMSD(nm)" << endl;
        for(int i=0; i<data[0].size(); ++i){
            f << i*dt << " ";
            for(int j=0; j<data.size(); ++j) f << data[j][i] << " ";
            f << endl;
        }
        f.close();
    }

private:
    vector<vector<float>> data;
    Selection fit_sel;
    vector<Selection> rms_sel;
};


CREATE_COMPILED_PLUGIN(rms)

