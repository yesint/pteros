#include "pteros/pteros_core.h"

using namespace std;
using namespace pteros;

int main(int argc, char** argv){
	try {
	    cout << "Pteros demo template program" << endl;
	    System sys("2lao.pdb");
	    cout << "Number of atoms is " << sys.num_atoms() << endl;
	} catch(Pteros_error e) {
	    e.print();
	}
}
