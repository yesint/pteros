#ifndef COMPILED_PLUGIN_H
#define COMPILED_PLUGIN_H

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
    class_<_name>("Task", init<pteros::Trajectory_processor*,pteros::Options_tree*>()) \
    .def_readwrite("label",&_name::label) \
    ; \
}

#endif
