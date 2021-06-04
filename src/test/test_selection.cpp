//auto INTEGER = g.rule() << lit('-')('?') << any_char("0123456789")('+');
//auto FLOAT = g.rule() << INTEGER << (any_char("eE") << INTEGER)('?')

#include <pteros/pteros.h>
#include <pteros/core/file_handler.h>
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


    string path="/home/semen/work/stored/Projects/ticagrelor/TIC_with_membranes/DOPC_TIC/with_NR12S";
    System s(path+"/confout.gro");
    //s.load(path+"/1500ns_wt_TIC.gro");
    s.load(path+"/traj_comp.xtc",0,10);

    Selection sel(s,"resname DOPC and within 0.3 pbc noself of resname TIC");
    sel.apply();
    cout << sel.coord_dependent() << " " << sel.size() <<  endl;

    for(int fr=0;fr<s.num_frames();++fr){
        sel.set_frame(fr);
        cout << sel.coord_dependent() << " " << sel.size() <<  endl;
    }

}
