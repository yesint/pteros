#include "pteros/pteros.h"
#include <Eigen/Core>

using namespace std;
using namespace pteros;
using namespace Eigen;

void help(){
    cout << "Usage:\n"
            "\t-solute <file>  - structure file with solute\n"
            "\t-solvent <file> - structure file with the periodic box of solvent\n"
            "\t\tDefaults to spc216.gro from Gromacs dir if Gromacs is installed\n"
            "\t\totherwise no default.\n"
            "\t-d <float>, default: 0.25 - minimal distance from solute to solvent in nm\n"
            "\t\tmeasured between the centers of atoms.\n"
            "\t-sel <string>, optional - custom selection of atoms to remove.\n"
            "\t\tExecuted after cutoff=d was applied.\n"
            "\t\tUseful for removing water from within lipid bilayer or protein cavities.\n"
            "\t\t**NOTE**: Always use 'pbc' in within selection to get meaningful result!\n"
            "\t\t**NOTE**: This selection is used 'as is', so be careful not to remove\n"
            "\t\tyour solute or doing other crazy things!\n"
            "\t-o <file>, default 'solvated.pdb' - output file."
         << endl;
}

int main(int argc, char* argv[]){
    try{
        cout << "===================================" << endl;
        cout << "==        pteros_solvent         ==" << endl;
        cout << "===================================" << endl;
        cout << "==  (C) Semen Yesylevskyy, 2016  ==" << endl;
        cout << "===================================" << endl;

        Options opt;
        parse_command_line(argc,argv,opt);

        if(opt.has("help")){
            help();
            return 0;
        }

        // Solute

        string solute_file = opt("solute").as_string();
        cout << "Loading solute from " << solute_file << " ..." << endl;
        System solute( solute_file );

        // Solvent

        System solvent;
        string solvent_file;
        // Look for $GMXDATA environmental variable
        if (const char* env_gmx = std::getenv("GMXDATA")) {
            solvent_file = opt("solvent",string(env_gmx)+"/top/spc216.gro").as_string();
        } else {
            solvent_file = opt("solvent").as_string();
        }
        cout << "Loading solvent from " << solvent_file << " ..." << endl;
        solvent.load( solvent_file );

        if(solvent.Box(0).is_triclinic())
            throw Pteros_error("Only rectangular solvent boxes are allowed!");

        // See how many solvent boxes should be used to cover solute

        // In arbitrary triclinic boxes find the maximal box coordinate
        Vector3f max_solute_coord = solute.Box(0).box_to_lab( solute.Box(0).extents() );
        Vector3f max_solvent_coord = solvent.Box(0).extents();

        Vector3i nbox;
        for(int i=0; i<3; ++i) nbox(i) = int(ceil(max_solute_coord(i)/max_solvent_coord(i)));

        cout << "Will use " << nbox.transpose() << " solvent boxes..." << endl;

        // Distribute solvent boxes
        {
            auto all = solvent.select_all();
            auto m = solvent.Box(0).get_matrix();
            solvent.distribute(all,nbox,m);
        }

        // Move min coords of solvent and solute to zero
        Vector3f solvent_min,solvent_max, solute_min, solute_max;
        auto solute_all = solute.select_all();
        auto solvent_all = solvent.select_all();
        solute_all.minmax(solute_min,solute_max);
        solvent_all.minmax(solvent_min,solvent_max);

        solvent_all.translate(-solvent_min);
        solute_all.translate(-solute_min);

        cout << "Finding solvent atoms outside the solute box..." << endl;

        // Cut solvent atoms outside the solute box
        vector<int> bad;        
        for(int i=0; i<solvent_all.size(); ++i){            
            if( !solute.Box(0).in_box(solvent_all.XYZ(i)) ) bad.push_back(solvent_all.Index(i));
        }

        // Select bad atoms
        Selection bad_sel(solvent,bad);
        // Select whole bad residues
        vector<Selection> bad_res;
        bad_sel.each_residue(bad_res);

        cout << "Found " << bad_res.size() << " solvent molecules outside the solute box..." << endl;
        for(auto& sel: bad_res){
            sel.set_beta(-1000);
        }

        // Find last index of solute
        int last_solute_ind = solute.num_atoms()-1;

        // append good solvent to solute
        solute.append(solvent("beta > -1000"));

        // select overlapping water
        float d = opt("d","0.25").as_float();
        string s = "by residue (index " + to_string(last_solute_ind+1)
                + "-" + to_string(solute.num_atoms()-1)
                + " and within "
                + to_string(d)+ " pbc of index 0-"
                + to_string(last_solute_ind) +")";

        Selection sel(solute, s);

        cout << "Found " << sel.size() << " overlaping solvent atoms at cutoff="
             << d <<"..."<< endl;

        // Remove overlapping water
        sel.set_beta(-1000);

        // If we have custom selection use it
        if(opt.has("sel")){            
            s = opt("sel").as_string();
            Selection sel(solute, s);
            cout << "Removing atoms from custom selection '" + s
                    + "' ("+to_string(sel.size())+" atoms)..." << endl;
            sel.set_beta(-1000);
        }

        // Translate back to initial box center
        solute_all.modify("beta > -1000");
        solute_all.translate(solute_min);

        // Report number of remaining solvent residues
        map<string,int> residues;
        int at=last_solute_ind+1;
        do {
            string resname = solute_all.Resname(at);
            int resind = solute_all.Resindex(at);

            // Find the end of this residue
            do {
                ++at;
            } while(solute_all.Resindex(at) == resind && at<solute_all.size());

            if(residues.count(resname)){
                // such resname is present
                ++residues[resname];
            } else {
                // new resname
                residues[resname] = 1;
            }

        } while(at<solute_all.size());

        cout << "Number of solvent molecules:" << endl;
        for(auto& it: residues){
            cout << it.first << " " << it.second << endl;
        }

        // Writing output
        solute_all.write( opt("o","solvated.pdb").as_string() );

    } catch(const Pteros_error& e) {
        e.print();
    }
}
