//auto INTEGER = g.rule() << lit('-')('?') << any_char("0123456789")('+');
//auto FLOAT = g.rule() << INTEGER << (any_char("eE") << INTEGER)('?')

#include <pteros/pteros.h>
#include <pteros/core/file_handler.h>
#include <pteros/core/grid.h>
#include <map>
#include <string>
#include <iostream>
#include <functional>
#include <memory>
#include <variant>
#include <charconv>
#include <string_view>

#include <unsupported/Eigen/CXX11/Tensor>

using namespace pteros;
using namespace std;
using namespace Eigen;


int main(int argc, char* argv[]){
    LOG()->set_level(spdlog::level::debug);

    Tensor<GridCell,3> t(3,5,6);
    t(1,1,1).clear();
    t(1,1,1).add_point(2,Vector3f{1,1,1});
    t(1,1,1).add_point(12,Vector3f{1,1,1});
    cout << t(1,1,1).size() << endl;

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
