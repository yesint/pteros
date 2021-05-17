#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "pteros/core/utilities.h"
#include "pteros/core/pteros_error.h"
#include "bindings_util.h"
#include "pteros/core/version.h"

namespace py = pybind11;
using namespace std;
using namespace pteros;
using namespace py::literals;

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
    m.doc() = "pteros bindings"; // module docstring

    make_bindings_Atom(m);
    make_bindings_System(m);
    make_bindings_Selection(m);
    make_bindings_Periodic_box(m);
    make_bindings_Frame(m);
    make_bindings_Distance_search(m);
    make_bindings_Options(m);
    make_bindings_Trajectory_reader(m);

    // Global stuff
    //m.attr("pteros_version") = py::cast(_git_revision);
    m.def("greeting",&greeting);

    py::class_<spdlog::logger,shared_ptr<spdlog::logger>>(m,"Logger")
        .def(py::init([](const string& name){            
            return create_logger(name);
        }))
        .def("info",[](spdlog::logger* log, const string& str){log->info(str);})
        .def("warn",[](spdlog::logger* log, const string& str){log->warn(str);})
        .def("error",[](spdlog::logger* log, const string& str){log->error(str);})
        .def("debug",[](spdlog::logger* log, const string& str){log->debug(str);})
        .def("critical",[](spdlog::logger* log, const string& str){log->critical(str);})
    ;
    m.def("set_log_level",&set_log_level);

    m.def("angle_between_vectors",&angle_between_vectors);
    m.def("project_vector",&project_vector);

    m.def("get_element_name",&get_element_name);
    m.def("get_element_number",&get_element_number);
    m.def("get_vdw_radius",&get_vdw_radius);

    m.def("rotation_transform",[](Vector3f_const_ref pivot, Vector3f_const_ref axis, float angle){
        return rotation_transform(pivot,axis,angle).matrix().transpose();
    });

    m.def("resname_1char",&resname_1char);
    m.def("resname_3char",&resname_3char);

    m.attr("noPBC") = py::cast(noPBC);
    m.attr("fullPBC") = py::cast(fullPBC);

    py::class_<Histogram>(m,"Histogram")
            .def(py::init<float,float,int>())
            .def("add",py::overload_cast<float,float>(&Histogram::add),"value"_a,"weight"_a=1.0)
            .def("add",py::overload_cast<const vector<float>&>(&Histogram::add))
            .def("add_cylindrical",&Histogram::add_cylindrical)
            .def("normalize",&Histogram::normalize,"norm"_a=0)
            .def("value",&Histogram::value)
            .def("position",&Histogram::position)
            .def("delta",&Histogram::delta)
            .def("get_bin",&Histogram::get_bin)
            .def_property_readonly("values",&Histogram::values)
            .def_property_readonly("positions",&Histogram::positions)
            .def_property_readonly("num_bins",&Histogram::num_bins)
            .def("save_to_file",&Histogram::save_to_file,"fname"_a,"x_shift"_a=0.0)
    ;

    py::class_<Histogram2D>(m,"Histogram2D")
            .def(py::init<float,float,int,float,float,int>())
            .def("add",&Histogram2D::add,"v1"_a,"v2"_a,"weight"_a=1.0)
            .def("normalize",&Histogram2D::normalize,"norm"_a=0)
            .def("save_to_file",&Histogram2D::save_to_file)
            .def("value",&Histogram2D::value)
    ;
}
