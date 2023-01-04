//auto INTEGER = g.rule() << lit('-')('?') << any_char("0123456789")('+');
//auto FLOAT = g.rule() << INTEGER << (any_char("eE") << INTEGER)('?')

#include <pteros/pteros.h>
#include "pteros/extras/voronoi_packing.h"
#include <ranges>

using namespace pteros;
using namespace std;
using namespace Eigen;


int main(int argc, char* argv[]){
    System s("/home/semen/work/siRNA/CG_Martini/SQ/CG_build/cgbuilder.gro");
}
