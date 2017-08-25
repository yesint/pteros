#include "pteros/core/selection.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

//PYBIND11_MAKE_OPAQUE(std::vector<int>);

namespace py = pybind11;
using namespace pybind11::literals;
using namespace std;
using namespace pteros;


PYBIND11_MODULE(pteros11, m) {
    m.doc() = "pybind11 pteros bindings"; // optional module docstring

    //py::bind_vector<std::vector<int>>(m, "VectorInt");
    using RowMatrixXf = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    py::class_<System>(m, "System")
        .def(py::init<const std::string &>())
        .def("num_atoms", &System::num_atoms)
        .def("num_frames", &System::num_frames)
        .def("getXYZ", py::overload_cast<int,int>(&System::XYZ, py::const_))
        .def("setXYZ", [](System& s,Vector3f_const_ref v,int i,int fr){ s.XYZ(i,fr)=v; })
        .def("load", py::overload_cast<string,int,int,int,std::function<bool(System*,int)>>(&System::load),
             "fname"_a, "b"_a=0, "e"_a=-1, "skip"_a=0, "on_frame"_a=nullptr)
        // Selecting
        .def("__call__", py::overload_cast<>(&System::operator()))
        .def("__call__", py::overload_cast<string,int>(&System::operator()),"str"_a,"fr"_a=0)
        .def("__call__", py::overload_cast<int,int>(&System::operator()))
        .def("__call__", py::overload_cast<const std::vector<int>&>(&System::operator()))
        .def("__call__", py::overload_cast<const std::function<void(const System&,int,std::vector<int>&)>&,int>(&System::operator()),"callback"_a,"fr"_a=0)
    ;

    py::class_<Selection>(m, "Selection")
        .def("size",&Selection::size)
        .def("__len__", &Selection::size)
        .def("get_index", &Selection::get_index)
        .def("translate", &Selection::translate)
        .def("get_xyz", [](Selection* sel){  RowMatrixXf m = sel->get_xyz().transpose(); return m; })
        .def("set_xyz", [](Selection* sel, MatrixXf_const_ref m){ sel->set_xyz(m.transpose()); })
    ;
}
