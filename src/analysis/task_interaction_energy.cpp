/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/

#include "task_interaction_energy.h"
#include <fstream>


using namespace std;
using namespace pteros;

Task_interaction_energy::Task_interaction_energy(Trajectory_processor* engine, Options_tree* opt):
    Task_base(engine,opt)
{   
    cutoff = options->get_value<float>("cut_off",0.25);
    is_periodic = options->get_value<bool>("periodic",false);
}

void Task_interaction_energy::pre_process(){
    data.clear();
    //searcher.create(cutoff,system,true,is_periodic);
}

bool Task_interaction_energy::process_frame(const Frame_info &info){    
    // See if we have self-energy
    (sel[0]==sel[1]) ? is_self_energy = true : is_self_energy = false;

    // Get contacts
    bon.clear();    
    if(!is_self_energy)
        Grid_searcher(cutoff,sel[0],sel[1],bon,true,is_periodic);
    else
        Grid_searcher(cutoff,sel[0],bon,true,is_periodic);
    // Get energy    
    Energy_components e = simulation->non_bond_energy(bon, sel[0].get_system()->Frame_data(sel[0].get_frame()));

/*
    cout << sel[0].get_text() << endl;
    cout << sel[1].get_text() << endl;
    Selection sel(system,"all");
    sel.set_frame(0);
    if(bon.size()){
        cout << info.absolute_frame << " " << bon.size() << endl;
        //cout << "***" << sel.XYZ(bon[0](0)).transpose() <<  " " << sel.XYZ(bon[0](1)).transpose() << endl;
        //cout << "***" << e.total << endl;
    }
*/

    data.push_back(e);
    mean.lj_14 += e.lj_14;
    mean.lj_sr += e.lj_sr;
    mean.q_14 += e.q_14;
    mean.q_sr += e.q_sr;
    mean.total += e.total;

    return true;
}

void Task_interaction_energy::post_process(const Frame_info& info){
    mean.lj_14 /= (float)info.valid_frame;
    mean.lj_sr /= (float)info.valid_frame;
    mean.q_14 /= (float)info.valid_frame;
    mean.q_sr /= (float)info.valid_frame;
    mean.total /= (float)info.valid_frame;

    // Output
    string fname = prefix+".dat";
    // Get time step in frames and time
    float dt = (info.last_time-info.first_time)/(float)(info.valid_frame);

    ofstream f(fname.c_str());
    f << "# Interaction energy of selections" << endl
      << "# '" << sel_name[0] << "' [" << sel_text[0] << "]" << endl
      << "# '" << sel_name[1] << "' [" << sel_text[1] << "]" << endl;
    f << "# Mean: " << mean.to_str() << endl;
    f << "# time total lj_sr lj_14 q_sr q_14:" << endl;
    for(int i=0; i<data.size(); ++i){
        f << i*dt << " " << data[i].to_str() << endl;
    }
    f.close();
}

//void Task_interaction_energy::window_started_slot(const Trajectory_processor_info& info){}
//void Task_interaction_energy::window_finished_slot(const Trajectory_processor_info& info){}

