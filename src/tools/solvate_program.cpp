#include "pteros/pteros.h"
#include <Eigen/Core>
#include "fmt/ostream.h"
#include "pteros/core/utilities.h"
#include "pteros/extras/solvate.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

string help(){
    return
R"(Usage:
-solute <file>  - structure file with solute
-solvent <file> - structure file with the box of solvent
    Defaults to spc216.gro from Gromacs dir if Gromacs is installed
    otherwise no default.
    Only rectangular solvent boxes are supported.
-d <float>, default: 0.25 - minimal distance from solute to solvent in nm
    measured between the centers of atoms.
-sel <string>, optional - custom selection of atoms to remove.
    Executed after cutoff=d was applied.
    Useful for removing water from within lipid bilayer or protein cavities.
    **NOTE**: Always use 'pbc' in within selection to get meaningful result!
    **NOTE**: This selection is used 'as is', so be careful not to remove
    your solute or doing other crazy things!
-o <file>, default 'solvated.pdb' - output file.
)";
}

int main(int argc, char* argv[]){
    try{
        greeting("pteros_solvate");

        LOG()->set_pattern("(%l)\t%v");

        Options opt;
        parse_command_line(argc,argv,opt);

        if(opt.has("help")){
            cout << help();
            return 0;
        }

        // Solute
        string solute_file = opt("solute").as_string();
        LOG()->info("Loading solute from '{}'...", solute_file);
        System solute( solute_file );

        string solvent_file;
        // Look for $GMXDATA environmental variable
        if (const char* env_gmx = std::getenv("GMXDATA")) {
            solvent_file = opt("solvent",string(env_gmx)+"/top/spc216.gro").as_string();
        } else {
            solvent_file = opt("solvent").as_string();
        }

        float d = opt("d","0.25").as_float();
        string custom_sel = opt("sel","").as_string();

        // Run solvate routine
        auto res = solvate(solute,d,solvent_file,custom_sel);

        // Writing output
        auto out = opt("o","solvated.pdb").as_string();
        LOG()->info("Writing output to '{}'...", out);
        solute().write(out);

    } catch(const PterosError& e) {
        LOG()->error(e.what());
    }
}
