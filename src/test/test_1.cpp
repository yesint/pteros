//-----------------------------
#include <pybind11/pybind11.h>

using namespace std;
namespace py = pybind11;

class A {
public:
   A(){}
   A(const A& other){}
   virtual void func()=0;
   virtual A* clone() const =0;
};

class A_trampoline: public A {
    using A::A;

    void func() override {
        PYBIND11_OVERLOAD_PURE(
            void, /* Return type */
            A,      /* Parent class */
            func,          /* Name of function in C++ (must match Python name) */
                          /* Argument(s) */
        );
    }

    A* clone() const override {
        PYBIND11_OVERLOAD_PURE(
            A*, /* Return type */
            A,      /* Parent class */
            clone,          /* Name of function in C++ (must match Python name) */
                          /* Argument(s) */
        );
    }

    //virtual A_trampoline* clone() const {
    //    return new A_trampoline(*this);
    //}

    //virtual A_trampoline* clone() const {
    //    auto this_object = py::cast(this);
    //    return py::cast<A_trampoline*>(this_object.get_type()(this_object).release());
    //}
};

// function dependent on clone functionality
void work_on_many_instances(const shared_ptr<A>& inst){
    std::vector<shared_ptr<A>> tasks;
    for(int i=0; i<100; ++i) tasks.push_back( shared_ptr<A>(inst->clone()) );
    for(int i=0; i<100; ++i) tasks[i]->func(); // Call python overload for all instances
}


PYBIND11_MODULE(module1, m) {

    py::class_<A,A_trampoline,std::shared_ptr<A>>(m, "A")
            .def(py::init_alias<>())
            .def(py::init_alias<const A_trampoline &>())
            .def("func",&A::func)
            .def("clone",&A::clone,py::return_value_policy::take_ownership)
    ;

    // Accepts PyObject with class derived in python from A
    m.def("work_on_many_instances", [](const py::object& o){
            auto p = o.cast<shared_ptr<A>>();
            //p->func();
            work_on_many_instances(p); // UPS! clone will try to create A, not the derived class...
       });
}
