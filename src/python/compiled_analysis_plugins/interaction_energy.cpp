/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2014, Semen Yesylevskyy
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
#include "pteros/core/grid_search.h"

using namespace std;
using namespace pteros;

class interaction_energy: public pteros::Compiled_plugin_base {
public:
    interaction_energy(pteros::Trajectory_processor* pr, pteros::Options_tree* opt): Compiled_plugin_base(pr,opt) {}

protected:

    void pre_process(){
        data.clear();
        cutoff = options->get_value<float>("cut_off",0.25);
        is_periodic = options->get_value<bool>("periodic",false);
        // Get selections
        std::list<string> sels = options->get_values<string>("selections");
        if(sels.size()<1 || sels.size()>2) throw Pteros_error("Either 1 or 2 selections should be passed");
        if(sels.size()==1){
            sel1.modify(system,  sels.front());
            sel2.modify(system,  sels.front());
        } else {
            sel1.modify(system,  sels.front());
            std::list<string>::iterator it = sels.begin();
            it++;
            sel2.modify(system, *it);
        }
    }

    void process_frame(const Frame_info &info){

        // See if we have self-energy
        (sel1==sel2) ? is_self_energy = true : is_self_energy = false;

        // Get contacts
        /*
        bon.clear();
        if(!is_self_energy)
            Grid_searcher(cutoff,sel1,sel2,bon,true,is_periodic);
        else
            Grid_searcher(cutoff,sel1,bon,true,is_periodic);             

        Energy_components e = system.non_bond_energy(bon, sel1.get_frame());
        */
        Energy_components e;
        for(int i=0;i<sel1.size(); ++i){
            for(int j=0;j<sel2.size(); ++j){
                system.add_non_bond_energy(e,sel1.Index(i),sel1.Index(j),0,true);
            }
        }

        data.push_back(e);
        mean.lj_14 += e.lj_14;
        mean.lj_sr += e.lj_sr;
        mean.q_14 += e.q_14;
        mean.q_sr += e.q_sr;
        mean.total += e.total;
    }

    void post_process(const Frame_info& info){
        mean.lj_14 /= (float)info.valid_frame;
        mean.lj_sr /= (float)info.valid_frame;
        mean.q_14 /= (float)info.valid_frame;
        mean.q_sr /= (float)info.valid_frame;
        mean.total /= (float)info.valid_frame;

        // Output
        string fname = label+".dat";
        // Get time step in frames and time
        float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);

        ofstream f(fname.c_str());
        f << "# Interaction energy of selections" << endl
          << "# [" << sel1.get_text() << "]" << endl
          << "# [" << sel2.get_text() << "]" << endl;
        f << "# Mean: " << mean.to_str() << endl;
        f << "# time total lj_sr lj_14 q_sr q_14:" << endl;
        for(int i=0; i<data.size(); ++i){
            f << i*dt << " " << data[i].to_str() << endl;
        }
        f.close();
    }

private:
    Selection sel1, sel2;
    bool is_self_energy;
    vector<Energy_components> data;
    Energy_components mean;
    vector<Eigen::Vector2i> bon;
    float cutoff;
    bool is_periodic;
};

CREATE_COMPILED_PLUGIN(interaction_energy)
