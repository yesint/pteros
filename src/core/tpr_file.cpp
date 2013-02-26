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

#include "tpr_file.h"
#include "pteros/core/pteros_error.h"
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

using namespace std;
using namespace pteros;
using namespace Eigen;


TPR_file::TPR_file(string fname, char mode): Mol_file(fname, mode)
{
    dump_name = "_dump_of_tpr_for_pteros";
    open_mode = mode;
    file_name = fname;
}

TPR_file::~TPR_file(){
    remove(dump_name.c_str());
}

bool TPR_file::do_read(System *sys, Frame *frame, Mol_file_content what){
    if(open_mode=='r'){

        cout << "Running gmxdump to dump content of TPR file " << file_name << "..." << endl;
        string cmd("gmxdump -sys -s ");
        cmd += file_name + " > " + dump_name;
        int ret = system(cmd.c_str());
        if(ret>0) throw Pteros_error("Error executing gmxdump! Check output above for error message.");
    } else {
        throw Pteros_error("TPR files could not be written from Pteros!");
    }

    // Start processing the dump
    cout << "Processing the dump..." << endl;

    ifstream f(dump_name.c_str());
    string line, dum, mode;
    stringstream ss;
    int natoms;
    boost::regex re;
    boost::cmatch m;

    while(getline(f,line)){
        if(line.size()<=1) continue; //Skip empty lines
        ss.clear(); ss.str(line);
        ss >> dum;

        if(dum=="atoms:"){
            getline(f,line);

            re.assign("(\\d+)");
            regex_search(line.c_str(), m, re);

            natoms = boost::lexical_cast<int>( string(m[1].first,m[1].second) );
            cout << "There are " << natoms << " atoms in topology" << endl;

            int max_type = 0; // The maximal atom type seen
            // Reading atoms
            cout << "Reading m,q,resid,type ..." << endl;
            re.assign("type=\\s*(\\d+),.+m=\\s*(\\S+),.+q=\\s*(\\S+),.+resind=\\s*(\\d+),");
            for(int i=0;i<natoms;++i){
                Atom at;

                getline(f,line);
                regex_search(line.c_str(), m, re);

                at.type = boost::lexical_cast<int>( string(m[1].first,m[1].second) );
                at.mass = boost::lexical_cast<float>( string(m[2].first,m[2].second) );
                at.charge = boost::lexical_cast<float>( string(m[3].first,m[3].second) );
                at.resid = boost::lexical_cast<int>( string(m[4].first,m[4].second) ) + 1;

                if(at.type>max_type) max_type = at.type;

                // Add new atom to the system
                append_atom_in_system(*sys,at);
            }

            // Now we know how many distinct atom types we have: max_type+1
            ff_in_system(*sys).nb_interactions.resize(max_type+1,max_type+1);

            getline(f,line);
            // Read atom names
            cout << "Reading names..." << endl;
            re.assign("\\\"(\\S+)\\\"");
            for(int i=0;i<natoms;++i){
                getline(f,line);
                regex_search(line.c_str(), m, re);
                atom_in_system(*sys,i).name = string(m[1].first,m[1].second);
            }

            getline(f,line);
            // Read atom typenames
            cout << "Reading type_names..." << endl;
            re.assign("name=\\\"(\\S+)\\\",");
            for(int i=0;i<natoms;++i){
                getline(f,line);
                regex_search(line.c_str(), m, re);
                atom_in_system(*sys,i).type_name = string(m[1].first,m[1].second);
            }

            // Reading number of residues
            getline(f,line);
            re.assign("(\\d+)");
            regex_search(line.c_str(), m, re);
            int Nres = boost::lexical_cast<int>( string(m[1].first,m[1].second) );
            cout << "There are " << Nres << " residues" << endl;

            // Reading residues
            cout << "Reading residue names..." << endl;
            re.assign("name=\\\"(\\S+)\\\",.+nr=(\\d+)");
            int cur_ind = 0;
            for(int i=0;i<Nres;++i){
                getline(f,line);
                regex_search(line.c_str(), m, re);
                int r = boost::lexical_cast<int>( string(m[2].first,m[2].second) );
                while(atom_in_system(*sys,cur_ind).resid==i){
                    atom_in_system(*sys,cur_ind).resname = string(m[1].first,m[1].second);
                    ++cur_ind;
                }
            }


        } else if(dum=="cgs:"){
            getline(f,line);
            re.assign("(\\d+)");
            regex_search(line.c_str(), m, re);
            int Ncgs = boost::lexical_cast<int>( string(m[1].first,m[1].second) );

            ff_in_system(*sys).charge_groups.resize(Ncgs);

            cout << "Reading charge groups..." << endl;
            re.assign("(\\d+)\\.\\.(\\d+)");
            for(int i=0;i<Ncgs;++i){
                getline(f,line);
                regex_search(line.c_str(), m, re);
                ff_in_system(*sys).charge_groups[i](0) =
                        boost::lexical_cast<int>( string(m[1].first,m[1].second) );
                ff_in_system(*sys).charge_groups[i](1) =
                        boost::lexical_cast<int>( string(m[2].first,m[2].second) );
                //cout << ff_in_system(sys).charge_groups[i].transpose() << endl;
            }
        } else if(dum=="excls:"){
            cout << "Reading exclusions..." << endl;
            getline(f,line);
            re.assign("(\\d+)");
            regex_search(line.c_str(), m, re);
            int Nexcl = boost::lexical_cast<int>( string(m[1].first,m[1].second) );

            getline(f,line);

            re.assign("\\{(.+)"); // From { to the end of string
            for(int i=0;i<Nexcl;++i){
                getline(f,line);

                regex_search(line.c_str(), m, re);
                string s(string(m[1].first,m[1].second));                
                cout << s << endl;

                // Parse line
                boost::regex spl(",\\s*|\\}\\s*");
                boost::sregex_token_iterator ii(s.begin(), s.end(), spl, -1);
                boost::sregex_token_iterator jj;
                while(ii != jj){
                    cout << *ii++ << endl;
                }

            }
        }

    }
    f.close();
}
