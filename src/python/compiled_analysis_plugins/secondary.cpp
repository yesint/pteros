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

using namespace std;
using namespace pteros;

class secondary: public Compiled_plugin_base {
public:
    secondary(Trajectory_processor* pr, Options_tree* opt): Compiled_plugin_base(pr,opt) {}

    string help(){
        return  "Purpose:\n"
                "\tComputes DSSP secondary structure for the system\n"
                "Output:\n"
                "\tFile <label>.dat containing the following columns:\n"
                "\ttime,DSSP\n"
                "\tThere is no space after ','! Spaces are DSSP codes themselves.\n"
                "Options:\n"
                "\tNone";
    }

protected:
    void pre_process(){
        // Output
        string fname = label+".dat";
        f.open(fname.c_str());
        f << "#frame,DSSP_code_string. NO space after ','!" << endl;
    }

    void process_frame(const Frame_info &info){
        f << info.valid_frame << "," << system.dssp() << endl;
    }

    void post_process(const Frame_info &info){        
        f.close();
    }    

private:
    ofstream f;
};

CREATE_COMPILED_PLUGIN(secondary)
