#include "pteros/pteros.h"
#include "pteros/core/pteros_error.h"
#include <iostream>

using namespace std;
using namespace pteros;

int main(int argc, char** argv){
	try {
        LOG()->info("Pteros demo template program");
	    System sys("2lao.pdb");
        LOG()->info("Number of atoms is {}", sys.num_atoms());
	} catch(Pteros_error e) {
        LOG()->error(e.what());
	}
}
