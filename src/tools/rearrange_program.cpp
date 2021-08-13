#include "pteros/pteros.h"
#include <Eigen/Core>
#include "fmt/ostream.h"
#include "pteros/core/utilities.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

string help(){
    return
R"(Usage:
-begin <sel1 sel2 sel3...>  - begin with these selections in the order of apperance
    Selections should not overlap.
-end <sel1 sel2 sel3...>  - end with these selections in the order of apperance
    Selections should not overlap.

Either "-begin" or "-end" or both should be specified.
All atoms not covered by "-begin" and "-end" are written "as is" between them.

-prefix <text>, optional - prefix to add to each selection (i.e. resname)
    Defaults to empty string.
    No space is added after prefix! Add one explicitly is needed.

-f <filename>, required - file to read
-o <filename>, optional, default: rearranged.pdb - file to write.
)";
}

int main(int argc, char* argv[]){
    try{
        greeting("pteros_rearrange");

        LOG()->set_pattern("(%l)\t%v");

        Options opt;
        parse_command_line(argc,argv,opt);

        if(opt.has("help")){
            cout << help();
            return 0;
        }

        vector<string> begin_sels = opt("begin","").as_strings();
        vector<string> end_sels = opt("end","").as_strings();

        if(begin_sels.size()==1 && begin_sels[0]=="") begin_sels.clear();
        if(end_sels.size()==1 && end_sels[0]=="") end_sels.clear();

        if(begin_sels.empty() && end_sels.empty())
            throw PterosError("You must specify at least some selections to rearrange!s");

        string prefix = opt("prefix","").as_string();
        if(prefix!=""){
            for(auto& s: begin_sels) s = prefix+s;
            for(auto& s: end_sels) s = prefix+s;
        }

        string fname = opt("f").as_string();
        System sys(fname);
        sys.rearrange(begin_sels,end_sels);
        sys().write(opt("o","rearranged.pdb").as_string());


    } catch(const PterosError& e) {
        LOG()->error(e.what());
    }
}
