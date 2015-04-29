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
#include "pteros/core/distance_search.h"
#include "pteros/core/system.h"

using namespace std;
using namespace pteros;

class interaction_energy: public pteros::Compiled_plugin_base {
public:
    interaction_energy(pteros::Trajectory_processor* pr, const pteros::Options& opt): Compiled_plugin_base(pr,opt) {}

    string help(){
        return  "Purpose:\n"
                "\tComputes non-bond interaction energy for each frame.\n"
                "\tIf one selection is provided computes self-energy.\n"
                "\tIf two selections are provided computes their interaction energy.\n"
                "\tCoordinate-dependent selections are updated for each frame.\n"
                "Output:\n"
                "\tFile <label>.dat containing the following columns:\n"
                "\ttime total lj_sr lj_14 q_sr q_14\n"
                "Options:\n"
                "\t-cutoff <float>, default: 0.25\n"
                "\t\tCutoff for energy computation\n"
                "\t-selections <string> [<string>]\n"
                "\t\tSelection texts for one or two selections\n"                
                "\t-periodic <bool>, default: false\n"
                "\t\tUse periodicity?\n";
    }
protected:

    void pre_process(){        
        cutoff = options("cutoff","0.25").as_float();
        is_periodic = options("periodic","false").as_bool();
        eps = options("eps","1.0").as_float();

        // Get selections
        std::vector<string> sels = options("selections").as_strings();
        if(sels.size()<1 || sels.size()>2) throw Pteros_error("Either 1 or 2 selections should be passed");
        if(sels.size()==1){            
            sel1.modify(system,sels.front());
            is_self_energy = true;            
        } else {            
            sel1.modify(system,sels[0]);
            sel2.modify(system,sels[1]);
            is_self_energy = false;            
        }

        // Output
        string fname = label+".dat";
        out.open(fname.c_str());

        if(is_self_energy){
            out << "# Interaction self-energy of selection" << endl
              << "# [" << sel1.get_text() << "]" << endl;
        } else {
            out << "# Interaction energy of selections" << endl
              << "# [" << sel1.get_text() << "]" << endl
              << "# [" << sel2.get_text() << "]" << endl;
        }

        out << "# time total lj_sr lj_14 q_sr q_14:" << endl;
    }

    void process_frame(const Frame_info &info){
        Energy_components e;

        if(is_self_energy){
            sel1.apply();
            e = sel1.non_bond_energy(cutoff,is_periodic);
        } else {
            sel1.apply();
            sel2.apply();
            e = non_bond_energy(sel1,sel2,cutoff,0,is_periodic);
        }

        out << info.absolute_time << " " << e.to_str() << endl;
    }

    void post_process(const Frame_info& info){        
        out.close();
    }

private:
    Selection sel1, sel2;
    bool is_self_energy;
    float cutoff;
    float eps;
    bool is_periodic;
    ofstream out;
};

CREATE_COMPILED_PLUGIN(interaction_energy)
