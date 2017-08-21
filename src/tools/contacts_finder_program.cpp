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

#include "contacts_finder.h"
#include "pteros/analysis/trajectory_reader.h"
#include <fstream>

using namespace pteros;
using namespace std;

int main(int argc, char* argv[]){

    try{
        Options options;

        parse_command_line(argc,argv,options);

        cout << "Creating trajectory processor..." << endl;
        Trajectory_reader processor(options);

        // Show help if asked
        if(argc==1){
            cout << processor.help() << endl;
            Contacts_finder::print_help();
            return 0;
        }

        cout << "Creating contacts finder..." << endl;
        processor.add_task(new Contacts_finder(options));

        // Do computation
        processor.run();


    } catch(const Pteros_error& e) {
        cout << "ERROR:" << e.what() << endl;
    }

}

