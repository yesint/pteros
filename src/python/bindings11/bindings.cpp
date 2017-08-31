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
void make_bindings_Options(py::module&);
void make_bindings_Trajectory_reader(py::module&);


PYBIND11_MODULE(_pteros11, m) {
    m.doc() = "pybind11 pteros bindings"; // optional module docstring

    make_bindings_Atom(m);
    make_bindings_System(m);
    make_bindings_Selection(m);
    make_bindings_Periodic_box(m);
    make_bindings_Frame(m);
    make_bindings_Distance_search(m);
    make_bindings_Options(m);
    make_bindings_Trajectory_reader(m);

    // Globas stuff
    py::class_<spdlog::logger>(m,"Logger")
        .def("info",[](spdlog::logger* log, const string& str){log->info(str);})
        .def("warn",[](spdlog::logger* log, const string& str){log->warn(str);})
        .def("error",[](spdlog::logger* log, const string& str){log->error(str);})
        .def("debug",[](spdlog::logger* log, const string& str){log->debug(str);})
        .def("critical",[](spdlog::logger* log, const string& str){log->critical(str);})
    ;

}
