#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "bindings_util.h"

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


PYBIND11_MODULE(_pteros, m) {
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
        .def(py::init([](const string& name){
            auto p = std::unique_ptr<spdlog::logger>(new spdlog::logger(name,Log::instance().console_sink));
            p->set_pattern(Log::instance().generic_pattern);
            return p;
        }))
        .def("info",[](spdlog::logger* log, const string& str){log->info(str);})
        .def("warn",[](spdlog::logger* log, const string& str){log->warn(str);})
        .def("error",[](spdlog::logger* log, const string& str){log->error(str);})
        .def("debug",[](spdlog::logger* log, const string& str){log->debug(str);})
        .def("critical",[](spdlog::logger* log, const string& str){log->critical(str);})
    ;

}
