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

    for(int i=0;i<1;++i){
        string path="/home/semen/work/stored/Projects/SquaMem/curved_charmm_sym/0_sym_flip";
        System s(path+"/last.gro");
        Selection sel(s,"within 0.5 of resid 1");
        LOG()->info("size {}",sel.size());

        vector<Vector2i> pairs;
        vector<float> dist;
        search_contacts(0.5,s("resid 1-10"),s("resid 11-20"),pairs,dist,true,fullPBC);
        LOG()->info("size {}",pairs.size());
    }

}
