/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2020, Semen Yesylevskyy
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



#include "rmsf.h"
#include "pteros/analysis/options.h"
#include "pteros/core/pteros_error.h"

using namespace std;
using namespace pteros;

int main(int argc, char** argv){

    try {
        Options options;
        parse_command_line(argc,argv,options);

        Trajectory_processor processor(options);

        // Show help if asked
        if(argc==1){
            cout << processor.help() << endl;
            RMSF::print_help();
            return 0;
        }

        std::unique_ptr<RMSF> engine(new RMSF(processor,options));

/*
        Selection s(sys,"all");
        Selection r(sys,"resid 1-100");

        cout << "-----> " << s.Resindex(0) << endl;

        Eigen::Vector3f v(0.1,0,0);
        for(int i=0;i<100;++i){
            sys.frame_dup(sys.num_frames()-1);
            r.set_frame(sys.num_frames()-1);
            r.translate(v);
            sys.Time(sys.num_frames()-1) = i+1;
        }

        s.write("zero.xtc");
*/

        processor.run();

        engine->save_results();

    } catch(const Pteros_error& e){
        cerr << "ERROR! " << e.what() << endl;        
    }
}


