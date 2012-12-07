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


#include "task_python.h"

using namespace pteros;
using namespace std;

int main(int argc, char* argv[]){

    try{
        /*
          This program reads options from command line as usual, but
          all options could also be set in python script and they will be overriden
          by the command line options. Command-line options has priority!
        */
        Options_tree options;

        options.from_command_line(argc,argv);

        // Show help
        if(options.count_options("help")){
            Trajectory_processor::print_help();            
            return 0;
        }


        if(argc==2 && !options.count_options("script")){
            // Try to interpret single cmd options as script name
            options.add_value("script",string(argv[1]));
        }


        cout << "Creating trajectory processor..." << endl;       
        Trajectory_processor processor;
        cout << "Creating python task driver..." << endl;
        // This will also parse options from python script and put then into options
        Task_python task(&processor,&options);
        processor.set_options(options);                

        // Do computation
        processor.run();

    } catch(Pteros_error e) {
        e.print();
    }

}
