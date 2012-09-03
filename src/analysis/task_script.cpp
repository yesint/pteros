#include "task_script.h"
#include <fstream>

using namespace std;
using namespace pteros;

Task_script::Task_script(Trajectory_processor* engine, Options_tree* opt,
            boost::shared_ptr<chaiscript::ChaiScript>& interp): Task_base(engine,opt) {

    using namespace chaiscript;

    chai = interp;
    chai->eval_file( opt->get_value<string>("script_file") );
    task = chai->eval("Task()");

    getter_name = "get_task"+boost::lexical_cast<string>(id);

    chai->add(fun<Boxed_Value()>(boost::bind(&Task_script::get_task,this)),getter_name);
    //chai->eval("print("+getter_name+"().selections_required())");

    chai->eval(getter_name+"().id = "+boost::lexical_cast<string>(id));

    // Add Frame_info to interpreter
    chai->add(user_type<Frame_info>(),"Frame_info");
    chai->add(fun(&Frame_info::absolute_frame),"absolute_frame");
    chai->add(fun(&Frame_info::valid_frame),"valid_frame");
}

void Task_script::pre_process(){
    chai->add(chaiscript::var(system), "sys");
    chai->eval(getter_name+"().system := sys");
    chai->add(chaiscript::var(&prefix),"pref");
    chai->eval(getter_name+"().prefix := pref");

    chai->eval(getter_name+"().pre_process()");
}

bool Task_script::process_frame(const Frame_info &info){
    chai->add(chaiscript::var(&info),"inf");
    bool ok = chai->eval<bool>(getter_name+"().process_frame(inf)");
    //chai->eval("print(\""+boost::lexical_cast<string>(id)+"\")");
    return ok;
}

void Task_script::post_process(const Frame_info &info){
    chai->add(chaiscript::var(info),"inf");
    chai->eval(getter_name+"().post_process(inf)");
}

//void Task_box::window_started_slot(const Trajectory_processor_info& info){}
//void Task_box::window_finished_slot(const Trajectory_processor_info& info){}

void Task_script::print_help(){
    cout << "Task script:\n"
            "-----------------\n"
            "Executes ChaiScript scripts for every frame.\n"
         << endl;
}
