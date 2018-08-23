#include "pteros/pteros.h"
#include <Eigen/Core>
#include "spdlog/fmt/ostr.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

string help(){
    return
R"(Usage:
-sel <sel1 sel2 sel3...>  - text for selections to rearrange in order of apperance
    All atoms not covered by passed selections are written "as is" after selected atoms.
    Selections should not overlap.
-prefix <text> - prefix to add to each selection (i.e. resname)
    Defaults to empty string.
    No space is added after prefix! Add one explicitly is needed.
-f <filename> - file to read
-o <filename>, optional, default: rearranged.pdb - file to write.
)";
}

int main(int argc, char* argv[]){
    try{
        cout << "===================================" << endl;
        cout << "==        pteros_rearrange       ==" << endl;
        cout << "===================================" << endl;
        cout << "==  (C) Semen Yesylevskyy, 2018  ==" << endl;
        cout << "===================================" << endl;

        LOG()->set_pattern("(%l)\t%v");

        Options opt;
        parse_command_line(argc,argv,opt);

        if(opt.has("help")){
            cout << help();
            return 0;
        }

        vector<string> sels = opt("sel").as_strings();

        string prefix = opt("prefix","").as_string();
        for(auto& s: sels) s = prefix+s;

        string fname = opt("f").as_string();
        System sys(fname);
        sys.rearrange(sels);
        sys().write(opt("o","rearranged.pdb").as_string());


    } catch(const Pteros_error& e) {
        LOG()->error(e.what());
    }
}
