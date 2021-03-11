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
#include <map>
#include "pteros/core/distance_search.h"
#include "pteros/core/system.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


TASK_PARALLEL(energy_par)
public:

    string help() override {
        return
R"(Purpose:
    Computes non-bond interaction energy for each frame.
    If one selection is provided computes self-energy.
    If two selections are provided computes their interaction energy.
    Coordinate-dependent selections are updated for each frame.
Output:
    File ebergy_<id>.dat containing the following columns:
    time total q lj
Options:
    -cutoff <float>, default: value from force field
        Cutoff for energy computation
    -sel <string> [<string>]
        Selection texts for one or two selections
    -periodic <bool>, default: true
        Use periodicity?
)";
    }
protected:

    void before_spawn() override {
        if(!system.force_field_ready()) throw PterosError("Need valid force field to compute energy!");

        cutoff = options("cutoff","0").as_float();
        is_periodic = options("periodic","true").as_bool();

        // Get selection texts
        sel_texts = options("sel").as_strings();
        if(sel_texts.size()<1 || sel_texts.size()>2) throw PterosError("Either 1 or 2 selections should be passed");
        is_self_energy = (sel_texts.size()==1) ? true : false;
    }

    void pre_process() override {
        if(is_self_energy){
            sel1.modify(system,sel_texts[0]);
        } else {
            sel1.modify(system,sel_texts[0]);
            sel2.modify(system,sel_texts[1]);
            if(check_selection_overlap({sel1,sel2})) throw PterosError("Selections could not overlap!");
            log->debug(sel1.get_text());
        }
    }

    void process_frame(const FrameInfo &info) override {
        Vector2f e;

        if(is_self_energy){
            sel1.apply();
            e = sel1.non_bond_energy(cutoff,is_periodic);
        } else {
            sel1.apply();
            sel2.apply();
            e = non_bond_energy(sel1,sel2,cutoff,0,is_periodic);
        }

        data[info.absolute_time] = e;
    }

    void post_process(const FrameInfo& info) override {
    }


    void collect_data(const std::vector<std::shared_ptr<TaskBase>>& tasks, int n_frames) override {
        for(const auto& it: tasks){
            auto h = dynamic_cast<energy_par*>(it.get());
            data.insert(h->data.begin(),h->data.end());
        }

        // Output
        ofstream out(fmt::format("energy_{}.dat",get_id()));

        if(is_self_energy){
            out << "# Interaction self-energy of selection" << endl
              << "# '" << sel1.get_text() << "'" << endl;
        } else {
            out << "# Interaction energy of selections" << endl
              << "# '" << sel1.get_text() << "'" << endl
              << "# '" << sel2.get_text() << "'" << endl;
        }

        out << "# time total q lj" << endl;
        out << "# cutoff: " << cutoff << endl;

        for(const auto& it: data){
            out << it.first << " " << it.second.sum() << " " << it.second.transpose() << endl;
        }

        out.close();
    }

private:
    Selection sel1, sel2;
    bool is_self_energy;
    float cutoff;    
    bool is_periodic;
    map<float,Vector2f> data;
    std::vector<string> sel_texts;
};

CREATE_COMPILED_PLUGIN(energy_par)




