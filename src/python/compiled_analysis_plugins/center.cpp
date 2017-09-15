/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;

TASK_SERIAL(center)
public:

    string help() override {
        return  "Purpose:\n"
                "\tComputes center of mass or geometric center of selection for each frame.\n"
                "\tIf selection is coordinate-dependent updates it every frame.\n"
                "Output:\n"
                "\tFile <label>.dat containing the following columns:\n"
                "\ttime center_x center_y center_z\n"
                "Options:\n"
                "\t--selection <string>\n"
                "\t\tSelection text\n"
                "\t--mass_weighted <true|false>, default: false\n"
                "\t\tCompute center of mass (true) or geometric center (false)";
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

    void process_frame(const Frame_info &info) override {
        sel.apply();
        f << info.absolute_time << " " << sel.center(use_mass).transpose() << endl;
        log->info("{} {}",info.absolute_time, sel.center(use_mass).transpose());
    }

    void post_process(const Frame_info &info) override {
        f.close();
    }    

private:
    Selection sel;
    bool use_mass;
    ofstream f;
};

CREATE_COMPILED_PLUGIN(center)

