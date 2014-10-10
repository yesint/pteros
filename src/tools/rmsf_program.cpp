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
