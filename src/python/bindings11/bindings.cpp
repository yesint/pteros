#include "pteros/core/selection.h"
#include <pybind11/pybind11.h>

//PYBIND11_MAKE_OPAQUE(std::vector<int>);

namespace py = pybind11;
using namespace std;
using namespace pteros;

// Forward declarations of all individual bindings of classes
void make_bindings_Atom(py::module&);
void make_bindings_System(py::module&);
void make_bindings_Selection(py::module&);
void make_bindings_Periodic_box(py::module&);
void make_bindings_Frame(py::module&);
void make_bindings_Distance_search(py::module&);


PYBIND11_MODULE(pteros11, m) {
    m.doc() = "pybind11 pteros bindings"; // optional module docstring

    make_bindings_Atom(m);
    make_bindings_System(m);
    make_bindings_Selection(m);
    make_bindings_Periodic_box(m);
    make_bindings_Frame(m);
    make_bindings_Distance_search(m);

    //py::bind_vector<std::vector<int>>(m, "VectorInt");

}
