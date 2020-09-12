//auto INTEGER = g.rule() << lit('-')('?') << any_char("0123456789")('+');
//auto FLOAT = g.rule() << INTEGER << (any_char("eE") << INTEGER)('?')

#include <pteros/pteros.h>
#include <pteros/core/mol_file.h>
#include <map>
#include <string>
#include <iostream>
#include <functional>
#include <memory>
#include <variant>

using namespace pteros;
using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]){
    LOG()->set_level(spdlog::level::debug);
    string path="/home/semen/work/stored/Projects/SquaMem/new_bending/cyl_DOPC/R100/";
    std::unique_ptr<Mol_file> h = Mol_file::open(path+"nojump.xtc",'r');

    int fr;
    float t;

    h->seek_frame(19);
    h->tell_current_frame_and_time(fr,t);
    cout << fr << " " << t << endl;

    h->seek_time(128010);
    h->tell_current_frame_and_time(fr,t);
    cout << fr << " " << t << endl;

}
