#include "bindings_util.h"

namespace py = pybind11;

// Forward declarations of all individual bindings of classes
void make_bindings_Membrane(py::module&);
void make_bindings_solvate(py::module&);

#ifdef WITH_OPENBABEL
void make_bindings_substructure_search(py::module&);
#endif

PYBIND11_MODULE(_pteros_extras, m) {
    m.doc() = "pteros extras bindings"; // module docstring

    make_bindings_Membrane(m);
    make_bindings_solvate(m);

#ifdef WITH_OPENBABEL
    make_bindings_substructure_search(m);
#endif

}
