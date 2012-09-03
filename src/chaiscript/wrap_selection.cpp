#include <Eigen/Core>
#include <chaiscript/chaiscript.hpp>
#include "pteros/core/system.h"
#include "pteros/core/selection.h"

using namespace chaiscript;
using namespace std;
using namespace Eigen;
using namespace pteros;

//void Selection_assign(Selection& lhs, Selection& rhs){ lhs = rhs; }

void wrap_selection(boost::shared_ptr<ChaiScript>& chai){

    chai->add(bootstrap::standard_library::vector_type<std::vector<Selection> >("Vector_Selection"));
    chai->add(bootstrap::standard_library::vector_type<std::vector<char> >("Vector_char"));
    chai->add(bootstrap::standard_library::vector_type<std::vector<string> >("Vector_string"));

    chai->add(user_type<Selection>(), "Selection");
    chai->add(constructor<Selection()>(), "Selection");
    chai->add(constructor<Selection(System&,string)>(), "Selection");
    chai->add(constructor<Selection(System&)>(), "Selection");
    chai->add(constructor<Selection(System&,int,int)>(), "Selection");

    chai->add(constructor<Selection(const Selection&)>(), "Selection");
    chai->add(fun(&Selection::operator=), "=");
    chai->add(fun(&Selection::operator==), "==");
    chai->add(fun(&Selection::operator!=), "!=");
    chai->add(fun<void(Selection::*)(Selection&)>(&Selection::append), "append");
    chai->add(fun<void(Selection::*)(int)>(&Selection::append), "append");
    chai->add(fun<void(Selection::*)(System&,string)>(&Selection::modify), "modify");
    chai->add(fun<void(Selection::*)(string)>(&Selection::modify), "modify");
    chai->add(fun<void(Selection::*)(System&,int,int)>(&Selection::modify), "modify");
    chai->add(fun<void(Selection::*)(int,int)>(&Selection::modify), "modify");
    chai->add(fun(&Selection::apply), "apply");
    chai->add(fun(&Selection::update), "update");
    chai->add(fun(&Selection::get_frame), "get_frame");
    chai->add(fun(&Selection::set_frame), "set_frame");
    chai->add(fun(&Selection::size), "size");
    chai->add(fun(&Selection::num), "num");
    chai->add(fun(&Selection::clear), "clear");
    chai->add(fun(&Selection::each_residue), "each_residue");
    chai->add(fun(&Selection::get_system), "get_system");
    chai->add(fun(&Selection::get_text), "get_text");
    chai->add(fun(&Selection::get_index), "get_index");
    chai->add(fun(&Selection::get_chain), "get_chain");
    chai->add(fun<void(Selection::*)(const vector<char>&)>(&Selection::set_chain), "set_chain");
    chai->add(fun<void(Selection::*)(char)>(&Selection::set_chain), "set_chain");
    chai->add(fun(&Selection::get_unique_chain), "get_unique_chain");
    chai->add(fun(&Selection::get_resid), "get_resid");
    chai->add(fun(&Selection::get_unique_resid), "get_unique_resid");
    chai->add(fun<void(Selection::*)(const vector<int>&)>(&Selection::set_resid), "set_resid");
    chai->add(fun<void(Selection::*)(int)>(&Selection::set_resid), "set_resid");
    chai->add(fun(&Selection::get_resindex), "get_resindex");
    chai->add(fun(&Selection::get_unique_resindex), "get_unique_resindex");
    chai->add(fun(&Selection::get_name), "get_name");
    chai->add(fun<void(Selection::*)(const vector<string>&)>(&Selection::set_name), "set_name");
    chai->add(fun<void(Selection::*)(string&)>(&Selection::set_name), "set_name");
    chai->add(fun<MatrixXf(Selection::*)()const>(&Selection::get_xyz), "get_xyz");
    chai->add(fun<void(Selection::*)(MatrixXf&)const>(&Selection::get_xyz), "get_xyz");
    chai->add(fun(&Selection::set_xyz), "set_xyz");
    // Overloads for get average
    chai->add(fun<MatrixXf(Selection*)>(boost::bind(&Selection::get_average,_1,0,-1)), "get_average");
    chai->add(fun<MatrixXf(Selection*,int)>(boost::bind(&Selection::get_average,_1,_2,-1)), "get_average");
    chai->add(fun<MatrixXf(Selection*,int,int)>(boost::bind(&Selection::get_average,_1,_2,_3)), "get_average");

    chai->add(fun(&Selection::get_mass), "get_mass");
    chai->add(fun<void(Selection::*)(const VectorXf&)>(&Selection::set_mass), "set_mass");
    chai->add(fun<void(Selection::*)(float)>(&Selection::set_mass), "set_mass");
}
