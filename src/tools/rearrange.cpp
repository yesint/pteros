#include "pteros/pteros.h"
#include <Eigen/Core>
#include "spdlog/fmt/ostr.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

void help(){
    cout << "Usage:\n"
            "\t-sel <sel1 sel2 sel3...>  - text for selections to rearrange in order of apperance\n"
            "\t-prefix <text> - prefix to add to each selection (i.e. resname)\n"
            "\t\tDefaults to empty string\n"
            "\t-f <filename> - file to read\n"
            "\t-o <filename>, optional, default: overwrite input file - file to write.\n"
         << endl;
}

int main(int argc, char* argv[]){
    try{
        cout << "===================================" << endl;
        cout << "==        pteros_rearrange       ==" << endl;
        cout << "===================================" << endl;
        cout << "==  (C) Semen Yesylevskyy, 2017  ==" << endl;
        cout << "===================================" << endl;

        LOG()->set_pattern("(%l)\t%v");

        Options opt;
        parse_command_line(argc,argv,opt);

        if(opt.has("help")){
            help();
            return 0;
        }

        vector<string> sels = opt("sel").as_strings();

        string prefix = opt("prefix","").as_string();
        for(auto& s: sels) s = prefix+s;

        string fname = opt("f").as_string();
        System sys(fname);
        sys.rearrange(sels);
        sys().write(opt("o",fname).as_string());


    } catch(const Pteros_error& e) {
        LOG()->error(e.what());
    }
}
