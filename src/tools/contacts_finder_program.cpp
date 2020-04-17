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



