#include <boost/algorithm/string/trim.hpp>

#define _CRT_SECURE_NO_WARNINGS

#include "pteros/bindings/chaiscript_engine.h"
#include "wrap_eigen.h"
#include "wrap_system.h"
#include "wrap_selection.h"
#include "pteros/core/system.h"
#include "pteros/core/selection.h"


#ifdef READLINE_AVAILABLE
#include <readline/readline.h>
#include <readline/history.h>
#else
char* readline(const char* p)
{
    std::string retval;
    std::cout << p ;
    std::getline(std::cin, retval);
#ifdef BOOST_MSVC
    return std::cin.eof() ? NULL : _strdup(retval.c_str());
#else
    return std::cin.eof() ? NULL : strdup(retval.c_str());
#endif
}
void add_history(const char*){}
void using_history(){}
#endif

void help(int n) {
    if ( n >= 0 ) {
        std::cout << "ChaiScript evaluator.  To evaluate an expression, type it and press <enter>." << std::endl;
        std::cout << "Additionally, you can inspect the runtime system using:" << std::endl;
        std::cout << "  dump_system() - outputs all functions registered to the system" << std::endl;
        std::cout << "  dump_object(x) - dumps information about the given symbol" << std::endl;
    } else {
        std::cout << "usage : chai [option]+" << std::endl;
        std::cout << "option:"                << std::endl;
        std::cout << "   -h | --help"         << std::endl;
        std::cout << "   -i | --interactive"  << std::endl;
        std::cout << "   -c | --command cmd"  << std::endl;
        std::cout << "   -v | --version"      << std::endl;
        std::cout << "   -    --stdin"        << std::endl;
        std::cout << "   filepath"            << std::endl;
    }
}

void version(int){
    std::cout << "chai: compiled " << __TIME__ << " " << __DATE__ << std::endl;
}

bool throws_exception(const boost::function<void ()> &f)
{
    try {
        f();
    } catch (...) {
        return true;
    }

    return false;
}

chaiscript::exception::eval_error get_eval_error(const boost::function<void ()> &f)
{
    try {
        f();
    } catch (const chaiscript::exception::eval_error &e) {
        return e;
    }

    throw std::runtime_error("no exception throw");
}

std::string get_next_command() {
    std::string retval("quit");
    if ( ! std::cin.eof() ) {
        char *input_raw = readline("eval> ");
        if ( input_raw ) {
            add_history(input_raw);
            retval = boost::trim_copy_if(std::string(input_raw),boost::is_any_of(" \t"));
            ::free(input_raw);
        }
    }
    if(   retval == "quit"
          || retval == "exit"
          || retval == "help"
          || retval == "version")
    {
        retval += "(0)";
    }
    return retval;
}

// We have to wrap exit with our own because Clang has a hard time with
// function pointers to functions with special attributes (system exit being marked NORETURN)
void myexit(int return_val) {
    exit(return_val);
}

void interactive(chaiscript::ChaiScript& chai)
{
    using_history();

    for (;;) {
        std::string input = get_next_command();
        try {
            // evaluate input
            chaiscript::Boxed_Value val = chai.eval(input);

            //Then, we try to print the result of the evaluation to the user
            if (!val.get_type_info().bare_equal(chaiscript::user_type<void>())) {
                try {
                    std::cout << chai.eval<boost::function<std::string (const chaiscript::Boxed_Value &bv)> >("to_string")(val) << std::endl;
                }
                catch (...) {} //If we can't, do nothing
            }
        }
        catch (const chaiscript::exception::eval_error &ee) {
            std::cout << ee.what();
            if (ee.call_stack.size() > 0) {
                std::cout << "during evaluation at (" << ee.call_stack[0]->start.line << ", " << ee.call_stack[0]->start.column << ")";
            }
            std::cout << std::endl;
        }
        catch (const std::exception &e) {
            std::cout << e.what();
            std::cout << std::endl;
        }
    }
}

void Frame_assign(pteros::Frame& lhs, pteros::Frame& rhs){ lhs = rhs; }
void Atom_assign(pteros::Atom& lhs, pteros::Atom& rhs){ lhs = rhs; }

boost::shared_ptr<chaiscript::ChaiScript> create_chaiscript_engine(int argc, char *argv[]){
    using namespace chaiscript;
    using namespace pteros;

    std::vector<std::string> usepaths;
    std::vector<std::string> modulepaths;   

    // Disable deprecation warning for getenv call.
#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable : 4996)
#endif

    const char *usepath = getenv("CHAI_USE_PATH");
    const char *modulepath = getenv("CHAI_MODULE_PATH");

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

    usepaths.push_back("");
    if (usepath)
    {
        usepaths.push_back(usepath);
    }

    modulepaths.push_back("");
    if (modulepath)
    {
        modulepaths.push_back(modulepath);
    }

    // Create engine
    boost::shared_ptr<ChaiScript> chai(new ChaiScript(modulepaths,usepaths));

    chai->add(fun(&myexit), "exit");
    chai->add(fun(&myexit), "quit");
    chai->add(fun(&help), "help");
    chai->add(fun(&version), "version");
    chai->add(fun(&throws_exception), "throws_exception");
    chai->add(fun(&get_eval_error), "get_eval_error");

    //-----------------
    // Add Eigen wrappers
    //-----------------
    wrap_eigen(chai);

    //-----------------
    // Adding System
    //-----------------
    wrap_system(chai);

    //-----------------
    // Adding Selection
    //-----------------
    wrap_selection(chai);

    //-----------------
    // Adding Frame
    //-----------------
    chai->add(user_type<Frame>(), "Frame");
    chai->add(constructor<Frame()>(), "Frame");
    chai->add(constructor<Frame(const Frame&)>(), "Frame");
    chai->add(fun<void(Frame&, Frame&)>(&Frame_assign), "=");

    //-----------------
    // Adding Atom
    //-----------------
    chai->add(user_type<Atom>(), "Atom");
    chai->add(constructor<Atom()>(), "Atom");
    chai->add(constructor<Atom(const Atom&)>(), "Atom");
    chai->add(fun<void(Atom&, Atom&)>(&Atom_assign), "=");


    return chai;
}

int run_chaiscript_engine(int argc, char *argv[], boost::shared_ptr<chaiscript::ChaiScript>& chai){

    for (int i = 0; i < argc; ++i) {
        if ( i == 0 && argc > 1 ) {
            ++i;
        }

        std::string arg( i ? argv[i] : "--interactive" );

        enum { eInteractive
               , eCommand
               , eFile
             } mode = eCommand ;

        if  ( arg == "-c" || arg == "--command" ) {
            if ( (i+1) >= argc ) {
                std::cout << "insufficient input following " << arg << std::endl;
                return EXIT_FAILURE;
            } else {
                arg = argv[++i];
            }
        } else if ( arg == "-" || arg == "--stdin" ) {
            arg = "" ;
            std::string line;
            while ( std::getline(std::cin, line) ) {
                arg += line + '\n' ;
            }
        } else if ( arg == "-v" || arg == "--version" ) {
            arg = "version(0)" ;
        } else if ( arg == "-h" || arg == "--help" ) {
            arg = "help(-1)";
        } else if ( arg == "-i" || arg == "--interactive" ) {
            mode = eInteractive ;
        } else if ( arg.find('-') == 0 ) {
            std::cout << "unrecognised argument " << arg << std::endl;
            return EXIT_FAILURE;
        } else {
            mode = eFile;
        }

        chaiscript::Boxed_Value val ;
        try {
            switch ( mode ) {
            case eInteractive : interactive(*chai); break;
            case eCommand     : val = chai->eval(arg); break;
            case eFile        : val = chai->eval_file(arg); break;
            }
        }
        catch (const chaiscript::exception::eval_error &ee) {
            std::cout << ee.what();
            if (ee.call_stack.size() > 0) {
                std::cout << "during evaluation at (" << *(ee.call_stack[0]->filename) << " " << ee.call_stack[0]->start.line << ", " << ee.call_stack[0]->start.column << ")";
                for (size_t j = 1; j < ee.call_stack.size(); ++j) {
                    if (ee.call_stack[j]->identifier != chaiscript::AST_Node_Type::Block
                            && ee.call_stack[j]->identifier != chaiscript::AST_Node_Type::File)
                    {
                        std::cout << std::endl;
                        std::cout << "  from " << *(ee.call_stack[j]->filename) << " (" << ee.call_stack[j]->start.line << ", " << ee.call_stack[j]->start.column << ")";
                    }
                }
            }
            std::cout << std::endl;
            return EXIT_FAILURE;
        }
        catch (std::exception &e) {
            std::cout << e.what() << std::endl;
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
