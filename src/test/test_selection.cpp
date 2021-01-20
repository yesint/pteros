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


    string path="/home/semen/work/Projects/pteros/github/pteros/src/test";
    System s(path+"/cg.gro");
    Selection sel1(s,"resid 3 ");
    Selection sel2(s,"resname  W");



    vector<Vector2i> bon;
    vector<float> dist;
    //search_contacts(1.5,sel1,bon,dist,true,noPBC);
    search_contacts(1.5,sel1,sel2,bon,dist,true,noPBC);

    cout << bon.size() << endl;

    Selection selw(s,"within 0.5 self nopbc of resid 3");
    cout <<selw.size() << endl;
    //cout << selw << endl;
    //for(int i=0;i<bon.size();++i)
    //    cout << bon[i].transpose() << " " << dist[i] << endl;

}
