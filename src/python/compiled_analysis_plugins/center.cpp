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
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;

TASK_SERIAL(center)
public:

    string help() override {
        return
R"(Purpose:
    Computes center of mass or geometric center of selection for each frame.
    If selection is coordinate-dependent updates it every frame.
Output:
    File <label>.dat containing the following columns:
    time center_x center_y center_z
Options:
    -selection <string>
        Selection text
    -mass_weighted <true|false>, default: false
        Compute center of mass (true) or geometric center (false)
)";
    }

protected:
    void pre_process() override {
        use_mass = options("mass","false").as_bool();
        string fname = fmt::format("center_id{}.dat",get_id());
        f.open(fname.c_str());
        f << "# time center_x center_y center_z" << endl;
        string sel_text = options("sel").as_string();
        sel.modify(system,sel_text);        
    }

    void process_frame(const FrameInfo &info) override {
        sel.apply();
        f << info.absolute_time << " " << sel.center(use_mass).transpose() << endl;
        log->info("{} {}",info.absolute_time, sel.center(use_mass).transpose());
    }

    void post_process(const FrameInfo &info) override {
        f.close();
    }    

private:
    Selection sel;
    bool use_mass;
    ofstream f;
};

CREATE_COMPILED_PLUGIN(center)





