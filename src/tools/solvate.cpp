#include "pteros/pteros.h"
#include <Eigen/Core>

using namespace std;
using namespace pteros;
using namespace Eigen;

void help(){
    cout << "Usage:\n"
            "\t-solute <file>  - structure file with solute\n"
            "\t-solvent <file> - structure file with the periodic box of solvent\n"
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

        Options opt;
        parse_command_line(argc,argv,opt);

        if(opt.has("help")){
            help();
            return 0;
        }

        cout << "Loading solute..." << endl;
        System solute( opt("solute").as_string() );
        cout << "Loading solvent..." << endl;
        System solvent( opt("solvent").as_string() );

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
        solvent.select_all().distribute(nbox,max_solvent_coord);

        // Move minimal coord solvent box to minimal coord of solute
        Vector3f solvent_min,solvent_max, solute_min, solute_max;
        Selection solute_all(solute,"all");
        Selection solvent_all(solvent,"all");
        solute_all.minmax(solute_min,solute_max);
        solvent_all.minmax(solvent_min,solvent_max);

        solvent_all.translate(solute_min-solvent_min);

        cout << "Finding solvent atoms outside the solute box..." << endl;

        // Cut solvent atoms outside the solute box
        vector<int> bad;
        Vector3f v;
        for(int i=0; i<solvent_all.size(); ++i){
            v = solute.Box(0).lab_to_box( solvent_all.XYZ(i)-solute_min );
            if(   v(0)>solute.Box(0).extent(0)
               || v(1)>solute.Box(0).extent(1)
               || v(2)>solute.Box(0).extent(2)
               || v(0)<0 || v(1)<0 || v(2) <0
              ) bad.push_back(solvent_all.Index(i));
        }

        cout << "Finding solvent residues outside the solute box..." << endl;

        // Select by residue for bad atoms
        string s("by residue index ");
        for(int i=0; i<bad.size(); ++i) s += to_string(bad[i])+" ";
        Selection bad_sel(solvent,s);

        cout << "Removing " << bad_sel.size() << " atoms outside the solute box..." << endl;
        solvent.remove(bad_sel);

        //solvent.Box(0) = solute.Box(0);
        //solvent.select_all().write("box.pdb");

        // Find last index of solute
        int last_solute_ind = solute.num_atoms()-1;

        // append solvent to solute
        solute.append(solvent);

        // select bad water
        float d = opt("d","0.25").as_float();
        s = "by residue (index " + to_string(last_solute_ind+1)
                + "-" + to_string(solute.num_atoms()-1)
                + " and within "
                + to_string(d)+ " pbc of index 0-"
                + to_string(last_solute_ind) +")";

        Selection sel(solute, s);

        cout << "Removing " << sel.size() << " overlaping solvent atoms at cutoff="
             << d <<"..."<< endl;

        // Remove bad water
        solute.remove(sel);

        // If we have custom selection use it
        if(opt.has("sel")){
            s = opt("sel").as_string();
            sel.modify(solute, s);
            cout << "Removing custom selection '" + s
                    + "' ("+to_string(sel.size())+" atoms)..." << endl;
            solute.remove(sel);
        }

        // Writing output
        solute.select_all().write( opt("o","solvated.pdb").as_string() );

    } catch(const Pteros_error& e) {
        e.print();
    }
}
