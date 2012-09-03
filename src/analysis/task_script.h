#ifndef TASK_SCRIPT_H
#define TASK_SCRIPT_H

#include "pteros/analysis/selection_analysis.h"

namespace pteros {


struct Task_script: public Task_base {

    Task_script(Trajectory_processor* engine, Options_tree* opt,
                boost::shared_ptr<chaiscript::ChaiScript>& interp);

    string task_name(){ return "script"; }
    int selections_required(){
        return chai->eval<int>(getter_name+"().selections_required()");
    }

    virtual void pre_process();
    virtual bool process_frame(const Frame_info& info);
    virtual void post_process(const Frame_info& info);
    //virtual void window_started_slot(const Trajectory_processor_info& info);
    //virtual void window_finished_slot(const Trajectory_processor_info& info);
    static void print_help();

    boost::shared_ptr<chaiscript::ChaiScript> chai;
    chaiscript::Boxed_Value task;
    chaiscript::Boxed_Value& get_task(){ return task; }
    string getter_name;
};

}

#endif
