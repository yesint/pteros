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

TASK_SERIAL(rms)
public:    

    string help() override {
        return
R"(Purpose:
    Computes RMSD of each frame for given selection.
    Selection should be coordinate-independent.
Output:
    File rms_id<id>.dat containing the following columns:
    time RMSD
    Also reports mean RMSD in the file header.
Options:
    -fit_sel <string>
        Fitting selection text
    -rms_sel <string string...>
        RMSD selections texts
    -nojump <distance>. Default: 0
        Remove jumps of atoms over periodic box boundary.
        Atoms, which should not jump, are unwrapped with
        given distance on the first frame.
        Zero means find unwrap distance automatically.
        Distance -1 means no jump removal
        this is only meaningful is very special cases!
)";
    }

protected:

    void pre_process() override {
        data.clear();        

        // rms_selections
        vector<string> strs = options("rms_sel").as_strings();
        if(strs.size()==0) throw PterosError("At least one rms selection required!");
        // Create selections
        for(auto& s: strs) rms_sel.emplace_back(system,s);

        if(check_selection_overlap(rms_sel)) throw PterosError("Selections should not overlap!");

        // fit_sel is optional
        string fit_str = options("fit_sel","").as_string();
        if(fit_str!=""){
            fit_sel.modify(system, fit_str);
        } else {
            log->info("Using first rms selection for fitting");
            fit_sel = rms_sel[0];
        }

        if(fit_sel.size()<3){
            throw PterosError("Can't fit selection with less than 3 atoms!");
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

    void process_frame(const pteros::FrameInfo &info) override {
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

    void post_process(const pteros::FrameInfo &info) override {
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




