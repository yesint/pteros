#include "pteros/bindings/chaiscript_engine.h"

using namespace chaiscript;
using namespace std;

int main(int argc, char *argv[])
{
    // Get instance of engine
    boost::shared_ptr<ChaiScript> chai = create_chaiscript_engine(argc,argv);

    // Run engine
    run_chaiscript_engine(argc,argv,chai);
}
