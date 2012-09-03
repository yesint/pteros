/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
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

#include "pteros/analysis/rmsf.h"
#include "pteros/pteros_core.h"

using namespace std;
using namespace pteros;

int main(int argc, char** argv){

    try {
        Options_tree options;
        options.from_command_line(argc,argv);
        //json_spirit::mValue m = engine->options.to_json();
        //json_spirit::write_stream(m,cout,true);
/*
        options.from_command_line(
"time ~/work/Projects/pteros/svn/pteros_build/src/analysis/pteros_rmsf "
"--trajectory_group "
"--range frame_range 0 1000 --end-range "
//"--window time_window 200 --end-window "
"/home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_noPBC_1.xtc "
"--end-trajectory_group "
"--structure_file /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/ref.pdb "
"--selection "
    "all "
    "--name all "
"--end-selection "
"--log_interval 1000 "
"--async true "
"--do_rmsd true"
            );
*/
        // Show help
        if(options.count_options("help") || argc==1){
            Trajectory_processor::print_help();
            RMSF::print_help();
            return 0;
        }

        Trajectory_processor processor(options);
        boost::shared_ptr<RMSF> engine(new RMSF(processor,options));

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

    } catch(Pteros_error e){
        cerr << "ERROR! " << e.what() << endl;        
    }
}
