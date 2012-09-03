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

#include "pteros/analysis/selection_analysis.h"

using namespace pteros;
using namespace std;

int main(int argc, char* argv[]){

    try{
        Options_tree options;

        options.from_command_line(argc,argv);
        /*
        string s =
        "--trajectory [ "
        "  /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_noPBC_1.xtc "
        "  /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_pdb2gmx.gro "
        "  /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/protein_only.top "
        "  --range frame_range 0 1000 "
        "] "
        "--async true "
        "--log_interval 100 "
        "--dump_input dump.txt "
        "--selection [ \"resid 364 to 528\" --name C1 ] "
        "--task rmsd "
        "--task box  "
        "--task [ center --weighted true] "
        "--task [interaction_energy --cut_off 0.25 --periodic false]";
        options.from_command_line(s);
        */

        // Show help
        if(options.count_options("help") || argc==1){
            Trajectory_processor::print_help();
            Selection_analysis::print_help();
            return 0;
        }

        cout << "Creating trajectory processor..." << endl;
        Trajectory_processor processor(options);
        cout << "Creating analysis driver..." << endl;
        Selection_analysis driver(processor,options);

        // Do computation
        processor.run();

    } catch(Pteros_error e) {
        e.print();
    }

}
