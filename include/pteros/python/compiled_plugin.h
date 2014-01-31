#ifndef COMPILED_PLUGIN_H
#define COMPILED_PLUGIN_H

#ifndef STANDALONE_PLUGINS

// Make a Python extension module from this plugin

#include "pteros/python/bindings_util.h"
#include "pteros/python/compiled_plugin_base.h"

using namespace boost::python;

#define CREATE_COMPILED_PLUGIN(_name) \
void Pteros_error_translator(const pteros::Pteros_error& e) { \
  PyErr_SetString(PyExc_UserWarning, const_cast<pteros::Pteros_error&>(e).what().c_str()); \
} \
BOOST_PYTHON_MODULE(_name) \
{ \
    import_array(); \
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray"); \
    register_exception_translator<pteros::Pteros_error>(&Pteros_error_translator); \
    class_<_name,boost::noncopyable>("Task", init<pteros::Trajectory_processor*,pteros::Options_tree*>()) \
    .def_readwrite("label",&_name::label) \
    .def("help",&_name::help) \
    ; \
}

#else //STANDALONE_PLUGINS

// Make a stand-alone executable from this plugin
#include "pteros/core/pteros_error.h"
#include "pteros/python/compiled_plugin_base.h"

using namespace pteros;
using namespace std;

#define CREATE_COMPILED_PLUGIN(_name) \
int main(int argc, char** argv){\
    try {\
        Options_tree options;\
        options.from_command_line(argc,argv);\
        Trajectory_processor engine(options);\
        _name task(&engine,&options);\
        task.label = #_name;\
        engine.run();\
    } catch(const Pteros_error& e) {\
        cout << e.what() << endl;\
    }\
}

#endif //STANDALONE_PLUGINS

#endif //COMPILED_PLUGIN_H
