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
#include "pteros/core/distance_search.h"
#include "pteros/core/system.h"

using namespace std;
using namespace pteros;


TASK_SERIAL(energy)
public:

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

    void pre_process() override {
        if(!system.force_field_ready()) throw Pteros_error("Need valid force fieled to compute energy!");

        cutoff = options("cutoff","0").as_float();
        is_periodic = options("periodic","false").as_bool();        

        // Get selections
        std::vector<string> sels = options("sel").as_strings();
        if(sels.size()<1 || sels.size()>2) throw Pteros_error("Either 1 or 2 selections should be passed");
        if(sels.size()==1){            
            sel1.modify(system,sels[0]);
            is_self_energy = true;            
        } else {            
            sel1.modify(system,sels[0]);
            sel2.modify(system,sels[1]);
            is_self_energy = false;            
        }

        // Output        
        out.open(fmt::format("energy_{}.dat",get_id()));

        if(is_self_energy){
            out << "# Interaction self-energy of selection" << endl
              << "# '" << sel1.get_text() << "'" << endl;
        } else {
            out << "# Interaction energy of selections" << endl
              << "# '" << sel1.get_text() << "'" << endl
              << "# '" << sel2.get_text() << "'" << endl;
        }

        out << "# time total q lj" << endl;
    }

    void process_frame(const Frame_info &info) override {
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

    void post_process(const Frame_info& info) override {
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

CREATE_COMPILED_PLUGIN(energy)




