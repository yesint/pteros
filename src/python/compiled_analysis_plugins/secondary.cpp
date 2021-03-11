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

using namespace std;
using namespace pteros;

TASK_SERIAL(secondary)
public:    

    string help() override {
        return
R"(Purpose:
    Computes DSSP secondary structure for the system
Output:
    File <label>.dat containing the following columns:
    time,DSSP
    There is no space after ','! Spaces are DSSP codes themselves.

    The file dssp_map_<label>.dat is then written with -map option.
Options:
    -onlynum <boolean>
        Output only the number of structured residues
    -map <boolean>
        Output secondary structure codes as numbers.
        All helices are encodes as 1, all beta sheets as 2 and all
        unstructured residues as 0.
        (useful for plotting)
)";
    }

protected:
    void pre_process() override {
        sel.modify(system, options("sel","protein").as_string());

        jump_remover.add_atoms(sel);


        onlynum = options("onlynum","false").as_bool();
        do_map = options("map","false").as_bool();

        // Output        
        f.open(options("out",fmt::format("dssp_{}.dat",get_id())).as_string());
        f << "#frame N :DSSP_code_string. NO space after ':'!" << endl;               
    }

    void process_frame(const FrameInfo &info) override {
        string s = sel.dssp();
        // Count all structured residues
        int N = std::count_if(s.begin(), s.end(), [](char c){return c!='T' && c!='S' && c!=' ';});
        f << info.valid_frame << " " << N;
        if(onlynum){
            f << endl;
        } else {
            f << " :" << s << endl;
        }

        if(do_map){
            // Convert all helices to 1, all beta sheets to 2, all the rest to 0
            vector<int> cur;
            for(char c: s){
                switch(c){
                case 'G':
                case 'H':
                case 'I':
                    cur.push_back(1);
                    break;
                case 'E':
                case 'B':
                    cur.push_back(2);
                    break;
                default:
                    cur.push_back(0);
                }
            }
            dssp_map.push_back(cur);
        }
    }

    void post_process(const FrameInfo &info) override {
        f.close();

        // Write map if asked
        if(do_map){
            f.open(options("out_map",fmt::format("dssp_map_{}.dat",get_id())).as_string());
            for(int fr=0;fr<dssp_map.size();++fr){
                for(int i=0;i<dssp_map[fr].size();++i) f << dssp_map[fr][i] << " ";
                f << endl;
            }
            f.close();
        }
    }    

private:
    bool onlynum;
    bool do_map;
    ofstream f;
    Selection sel;

    vector<vector<int>> dssp_map;
};

CREATE_COMPILED_PLUGIN(secondary)




