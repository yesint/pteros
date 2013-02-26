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

#include "pteros/analysis/selection_analysis.h"
#include "task_center.h"
#include "task_box.h"
#include "task_rmsd.h"
#include "task_interaction_energy.h"
#include "task_distance_matr.h"
#include "task_align_parts.h"
#include "task_covar_matr.h"

using namespace pteros;
using namespace Eigen;

void Task_base::setup(){
    // Create selections
    for(int i=0;i<sel_text.size();++i){
        sel.push_back(Selection());
        sel.back().modify(system,sel_text[i]);
        //cout << "\tmaking sel " << sel_text[i] << endl;
    }
    // Get file prefixes for output files
    prefix = "task."+boost::lexical_cast<string>(id)+"."+task_name();
    if(sel.size()){
        prefix += "-sel";
        for(int i=0; i<sel.size(); ++i){
            prefix += "."+sel_name[i];
        }
    }

}

Task_base::Task_base(Trajectory_processor *engine, Options_tree* opt):
    Consumer(engine)
{    
    options = opt;
    is_ready = false;
}

void Task_base::create(std::vector<std::string> &sel_texts, std::vector<std::string> &sel_names){
    if(sel_texts.size() != sel_names.size())
        throw Pteros_error("Mismatch beetwen names and selection texts!");
    if(selections_required()>=0 && sel_texts.size() != selections_required())
        throw Pteros_error("Wrong number of selections given!");

    sel_text = sel_texts;
    sel_name = sel_names;

    is_ready = true;
}

void Task_base::before_each_frame(){
    // Here we update all coordinate-dependent selections
    for(int i=0; i<sel.size(); ++i){
        sel[i].set_frame(0);
        // This will call apply() automatically for coordinate-dependent selections
    }
}


Selection_analysis::Selection_analysis(){
}


void Selection_analysis::create(Trajectory_processor& proc, Options_tree& opt){
    options = &opt;
    engine = &proc;

    vector<string> sel_names, sel_texts;        

    // Extract all tasks
    for(Options_tree* t: options->get_options("task")){
        // Determine type of task        
        string task_name = t->get_value<string>("");

        // Create a task using the factory
        tasks.push_back( task_factory(task_name,t) );

        sel_names.clear();
        sel_texts.clear();

        if(tasks.back()->selections_required() == 0){
            // Zero-selection tasks are special
            tasks.back()->create(sel_texts,sel_names);
            cout << "\tAdded task '" << task_name << endl;
        } else {
            // cycle over all selections
            int k = 0;
            for(Options_tree* s: options->get_options("selection")){
                // Get name of selection
                sel_names.push_back( s->get_value<string>("name",boost::lexical_cast<string>(k)) );
                // Get test of selection
                sel_texts.push_back( s->get_value<string>("") );

                if(tasks.back()->selections_required() == sel_names.size()){
                    // We added enough selections
                    if(tasks.back()->task_ready()) tasks.push_back( task_factory(task_name,t) );
                    tasks.back()->create(sel_texts,sel_names);

                    cout << "\tAdded task '" << task_name << "' for selections";
                    for(int i=0; i<sel_names.size();++i) cout << " '" << sel_names[i] << "'";
                    cout << endl;

                    sel_names.clear();
                    sel_texts.clear();
                }
                ++k;
            } // End of cycle over selections
        }
        // Tasks consuming all selections are also special
        if(tasks.back()->selections_required() == -1){
            tasks.back()->create(sel_texts,sel_names);

            cout << "\tAdded task '" << task_name << "' for selections";
            for(int i=0; i<sel_names.size();++i) cout << " '" << sel_names[i] << "'";
            cout << endl;
        }
    }
}

void Selection_analysis::print_help(){
    cout << "Specific options:\n"
            "-----------------\n"
            "--selection sel_text [--name sel_name]:\n"
            "\tRequired, may be given several times. Selection for analysis.\n"
            "\tsel_name is used for output files. If not specified numbers are used.\n"

            "--task task_name [...task sub-options...]:\n"
            "\tRequired, may be given several times.\n"
            "\tTask to perform on specified selections.\n"
            "\tTasks may require 0, 1 or 2 selections, which are taken sequencially\n"
            "\tfrom the list of provided seletcions. If 5 selection are given\n"
            "\tthen the task operating on 2 selections will be called for\n"
            "\tpair 1:2 and 3:4.\n"

            "Tasks:\n"
            "------\n"
         << endl;
    // Help for individual tasks
    Task_box::print_help();
    //Task_center.print_help();
    //Task_rmsd.print_help();
    //Task_interaction_energy.print_help();
}

boost::shared_ptr<Task_base> Selection_analysis::task_factory(std::string task_name, Options_tree* opt){
    if(task_name=="center"){
        return boost::shared_ptr<Task_base>(new Task_center(engine,opt));
    } else if(task_name=="box"){
        return boost::shared_ptr<Task_base>(new Task_box(engine,opt));
    } else if(task_name=="interaction_energy"){
        return boost::shared_ptr<Task_base>(new Task_interaction_energy(engine,opt));
    } else if(task_name=="distance_matr"){
        return boost::shared_ptr<Task_base>(new Task_distance_matr(engine,opt));
    } else if(task_name=="align_parts"){
        return boost::shared_ptr<Task_base>(new Task_align_parts(engine,opt));
    } else if(task_name=="covar_matr"){
        return boost::shared_ptr<Task_base>(new Task_covar_matr(engine,opt));
    } else {
        throw Pteros_error("Task '"+task_name+"' is not recognized!");
    }

}
