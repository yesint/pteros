//auto INTEGER = g.rule() << lit('-')('?') << any_char("0123456789")('+');
//auto FLOAT = g.rule() << INTEGER << (any_char("eE") << INTEGER)('?')

#include <pteros/pteros.h>
#include "pteros/extras/voronoi_packing.h"

using namespace pteros;
using namespace std;
using namespace Eigen;


int main(int argc, char* argv[]){
    System s("/home/semen/work/CLCK/CLCKa/r1.gro");
    vector<Selection> sels;
    sels.emplace_back(s("protein"));
    sels.emplace_back(s("resname POPC"));
    sels.emplace_back(s("resname CHL1"));
    sels.emplace_back(s("resname TIP3"));
    sels.emplace_back(s("resname SOD"));
    sels.emplace_back(s("resname CLA"));
    
    Voronoi3D voro(sels);
    voro.compute();
}
