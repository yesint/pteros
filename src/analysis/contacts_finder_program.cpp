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

#include "pteros/analysis/contacts_finder.h"

using namespace pteros;
using namespace std;

int main(int argc, char* argv[]){

    try{
        Options_tree options;

        options.from_command_line(argc,argv);
/*
        options.from_command_line("--trajectory_group "
            "/home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_noPBC_1.xtc "
            "--range frame_range 0 100 --end-range "            
            " --structure_file /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_pdb2gmx.gro "
            " --end-trajectory_group "
            " --selections \"resid 364 to 528\" "
            " \"resid 1 to 342 or resid 529 to 869\" "
            " --topology_file /home/semen/work/Projects/kornelyuk/Sasha/dimer_md/protein_only.top "
            " --method cut_off 0.25 "
            " --async true "
            "--log_interval 100 ");
*/

        // Show help
        if(options.count_options("help") || argc==1){
            Trajectory_processor::print_help();
            Contacts_finder::print_help();
            return 0;
        }

        cout << "Creating trajectory processor..." << endl;
        Trajectory_processor processor(options);
        cout << "Creating contacts finder..." << endl;
        boost::shared_ptr<Contacts_finder> finder(new Contacts_finder(processor,options));

        // Do computation
        processor.run();

        cout << "Saving results..." << endl;
        {
            ofstream ff("results.json");
            finder->save(ff);
            ff.close();
        }

        cout << "Saving per frame stats..." << endl;
        {
            ofstream ff("data.dat");
            finder->per_frame_stats(ff);
            ff.close();
        }

    } catch(Pteros_error e) {
        e.print();        
    }

}

