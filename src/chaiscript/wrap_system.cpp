#include <Eigen/Core>
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include <boost/bind.hpp>
#include <chaiscript/chaiscript.hpp>

using namespace chaiscript;
using namespace std;
using namespace Eigen;
using namespace pteros;

void System_atoms_dup1(System* sys, const std::vector<int>& ind){
    sys->atoms_dup(ind,NULL);
}

void System_atoms_dup2(System* sys, const std::vector<int>& ind, Selection* res_sel){
    sys->atoms_dup(ind,res_sel);
}

float System_distance1(System* sys, int i, int j, int fr){
    return sys->distance(i,j,fr,false);
}

float System_distance2(System* sys, int i, int j, int fr, bool is_periodic){
    return sys->distance(i,j,fr,is_periodic);
}

float System_distance3(System* sys, const Eigen::Vector3f& p1, const Eigen::Vector3f& p2, int fr){
    return sys->distance(p1,p2,fr,false);
}


float System_distance4(System* sys, const Eigen::Vector3f& p1, const Eigen::Vector3f& p2,
                       int fr, bool is_periodic){
    return sys->distance(p1,p2,fr,is_periodic);
}


void wrap_system(boost::shared_ptr<ChaiScript>& chai){
    chai->add(bootstrap::standard_library::vector_type<std::vector<int> >("Vector_int"));
    chai->add(bootstrap::standard_library::vector_type<std::vector<Atom> >("Vector_Atom"));
    chai->add(bootstrap::standard_library::vector_type<std::vector<Vector3f> >("Vector_Vector3f"));    

    chai->add(user_type<System>(), "System");
    chai->add(constructor<System()>(), "System");
    chai->add(constructor<System(string)>(), "System");
    chai->add(constructor<System(const System &)>(), "System");
    chai->add(fun(&System::num_atoms), "num_atoms");
    chai->add(fun(&System::num_frames), "num_frames");
    // Overloads for load
    chai->add(fun<void(System*,string)>( boost::bind(&System::load,_1,_2,0,-1,0) ), "load");
    chai->add(fun<void(System*,string,int)>( boost::bind(&System::load,_1,_2,_3,-1,0) ), "load");
    chai->add(fun<void(System*,string,int,int)>( boost::bind(&System::load,_1,_2,_3,_4,0) ), "load");
    chai->add(fun<void(System*,string,int,int,int)>( boost::bind(&System::load,_1,_2,_3,_4,_5) ), "load");

    chai->add(fun(&System::frame_dup), "frame_dup");
    chai->add(fun(&System::set_frame), "set_frame");
    chai->add(fun(&System::frame_copy), "frame_copy");
    // Overloads for frame_delete
    chai->add(fun<void(System*)>( boost::bind(&System::frame_delete,_1,0,-1) ), "frame_delete");
    chai->add(fun<void(System*,int)>( boost::bind(&System::frame_delete,_1,_2,-1) ), "frame_delete");
    chai->add(fun<void(System*,int,int)>( boost::bind(&System::frame_delete,_1,_2,_3) ), "frame_delete");

    chai->add(fun(&System::Frame_data), "Frame_data");
    // Overloads of clear
    chai->add(fun<void(System*)>( boost::bind(&System::clear,_1,false) ), "clear");
    chai->add(fun<void(System*,bool)>( boost::bind(&System::clear,_1,_2) ), "clear");

    chai->add(fun(&System::update_selections), "update_selections");
    chai->add(fun(&System::Box), "Box");
    chai->add(fun(&System::Time), "Time");
    chai->add(fun(&System::XYZ), "XYZ");
    chai->add(fun(&System::is_box_triclinic), "is_box_triclinic");
    chai->add(fun(&System::get_box_vectors_angles), "get_box_vectors_angles");
    chai->add(fun(&System::frame_append), "frame_append");
    chai->add(fun(&System::assign_resindex), "assign_resindex");
    // Overloads for atoms_dup
    // boost::bind is stupid and can't bind NULL arguments, so we need a manual overload
    chai->add(fun(&System_atoms_dup1), "atoms_dup");
    chai->add(fun(&System_atoms_dup2), "atoms_dup");
    //chai->add(fun<void(System*,const vector<int>&,Selection*)>( boost::bind(&System::atoms_dup,_1,_2,_3) ), "atoms_dup");
    chai->add(fun(&System::atoms_add), "atoms_add");
    // For some reason atoms_delete is not wrapped without explicit signature
    chai->add(fun<void(System*, const std::vector<int>&)>(&System::atoms_delete), "atoms_delete");
    chai->add(fun(&System::append),"append");
    // Distance overloads. Since there are two overloads each with 5 args
    // trying to bind them with boost::bind is a nightmare (it can't distinguish them).
    // Thus simpler to make wrappers
    chai->add(fun(&System_distance1),"distance");
    chai->add(fun(&System_distance2),"distance");
    chai->add(fun(&System_distance3),"distance");
    chai->add(fun(&System_distance4),"distance");

    chai->add(fun(&System::wrap_to_box),"wrap_to_box");

}
