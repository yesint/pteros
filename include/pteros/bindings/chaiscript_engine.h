#ifndef ENGINE_H
#define ENGINE_H

#include <chaiscript/chaiscript.hpp>

boost::shared_ptr<chaiscript::ChaiScript> create_chaiscript_engine(int argc, char *argv[]);
int run_chaiscript_engine(int argc, char *argv[], boost::shared_ptr<chaiscript::ChaiScript>& chai);

#endif
