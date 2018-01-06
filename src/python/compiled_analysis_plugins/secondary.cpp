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

using namespace std;
using namespace pteros;

class secondary: public Compiled_plugin_base {
public:
    secondary(Trajectory_processor* pr, const Options& opt): Compiled_plugin_base(pr,opt) {}

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
        f << info.valid_frame << "," << system.dssp(0) << endl;
    }

    void post_process(const Frame_info &info){        
        f.close();
    }    

private:
    ofstream f;
};

CREATE_COMPILED_PLUGIN(secondary)

