#include "bindings_util.h"

namespace py = pybind11;

// Forward declarations of all individual bindings of classes
void make_bindings_Membrane(py::module&);


PYBIND11_MODULE(_pteros_extras, m) {
    m.doc() = "pteros extras bindings"; // module docstring

    make_bindings_Membrane(m);
}
