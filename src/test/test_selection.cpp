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
#include <charconv>
#include <string_view>

using namespace pteros;
using namespace std;
using namespace Eigen;


int main(int argc, char* argv[]){
    LOG()->set_level(spdlog::level::debug);


    string path="/home/semen/work/stored/Projects/ticagrelor/TIC_with_membranes/realistic_membr_with_TIC_Florentin/data";
    System s(path+"/1500ns_wt_TIC.gro");

    s.write("a.xtc");

}
