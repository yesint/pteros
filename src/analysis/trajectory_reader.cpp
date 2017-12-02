/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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

#include "pteros/analysis/trajectory_reader.h"
#include "message_channel.h"
#include "data_container.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/mol_file.h"
#include "task_driver.h"
#include "traj_file_reader.h"
#include "pteros/core/logging.h"
#include <thread>

using namespace pteros;
using namespace std;


string Trajectory_reader::help(){
    return  "Trajectory processing options:\n"
            "General usage:\n"
            "\t-f filename1 filename2 ... <processing options>\n"
            "Files:\n"
            "\t* Exactly one structure file (PDB or GRO)\n"
            "\t  If not specified, topology PTTOP file must be given instead.\n"
            "\t* Topology PTTOP file (converted from Gromacs .tpr by tpr2pteros.py)\n"
            "\t  If structure file is also present only topology is read from this file.\n"
            "\t  If structure file is not present the coordinates are also read.\n"
            "\t* One or more trajectory files (TRR, XTC, TNG or DCD).\n"
            "\t  TNG files also contain the structure, so if no structure file\n"
            "\t  is given the structure is read from the first TNG file.\n"
            "\n"
            "\tFiles may appear in any order, but trajectory files will be processed\n"
            "\tin the order of their appearance.\n\n"

            "Processing options:\n"

            "\t-path <string>\n"
            "\t\toptional path which will be prepended to all data files, default: empty string\n"

            "\t-b <value[suffix]>\n"
            "\t\tbeginning of processing (starting frame or time), default: 0\n"

            "\t-e <value[suffix]>\n"
            "\t\tend of processing (end frame or time), default: -1 (up to the end)\n"

            "\t-skip <n>\n"
            "\t\tProcess only each n'th frame, default: -1 (process each frame)\n"

            "\t-t0 <t[suffix]>\n"
            "\t\tCustom starting time, default: -1 (use value from first frame)\n"
            "\t\tUseful if trajectory does not contain time stamps\n"
            "\t\tor if the starting time is incorrect.\n"
            "\t\tIf set and dt is not given sets dt to 1.0!\n"

            "\t-dt <t[suffix]>\n"
            "\t\tCutom time step, default: -1 (use value from trajectory)\n"
            "\t\tUseful if trajectory does not contain time stamps.\n"
            "\t\tIf set and start is not given sets start to 0.0!\n"

            "\t-log <n>\n"
            "\t\tPrints logging information on each n-th frame, default: -1 (no logging)\n"

            "\t-buffer <n>\n"
            "\t\tNumber of frames, which are kept in memory, default: 10\n"
            "\t\tOnly touch this if individual frames are very large.\n"

            "Suffixes:\n"
            "\tAll parameters marked as <value[suffix]> accept the following optional suffixes:\n"
            "\t\t(no suffix) - value is in frames\n"
            "\t\tfr - value is in frames\n"
            "\t\tt - value is time in picoseconds (value used as is)\n"
            "\t\tps - value is time in picoseconds (value used as is)\n"
            "\t\tns - value is time in nanoseconds (value multiplied by 10^3)\n"
            "\t\tus - value is time in microseconds (value multiplied by 10^6)\n"
            "\t\tms - value is time in milliseconds (value multiplied by 10^9)\n"
            "\tParameters marked as <t[suffix]> does not accept fr suffix.\n"
            "\tIn this case no suffix means ps.\n"
            ;
}



Trajectory_reader::Trajectory_reader()
{

}

Trajectory_reader::Trajectory_reader(const Options &opt): options(opt)
{

}

void Trajectory_reader::run(){    
    // Separate logger (not registered since only used here)
    auto log = create_logger("trj_processor");

    // Preparation stage
    log->debug("Starting trajectory processing");
    auto start = chrono::steady_clock::now();

    // Get file path
    auto file_path = options("path","").as_string();
    if(file_path!="" && file_path.back()!='/') file_path += "/";

    // Get files to work with
    auto file_list = options("f").as_strings();

    // Prepend path to all files
    for(auto& f: file_list) f = file_path + f;

    // Get all string parameters and find out how many different files we have
    // and arrange them structure->topology->traj
    string top_file = "";
    string structure_file = "";    

    for(string& s: file_list){                
        auto h = Mol_file::recognize(s);
        auto c = h->get_content_type();

        // traj file is always added as traj even if this is TNG
        if(c.traj()){
            traj_files.push_back(s);
            continue; // Avoid adding TNG twice also as structure file
        }

        if(c.top()){
            if(top_file=="")
                top_file = s;
            else
                throw Pteros_error("Only one topology file allowed!");
        }

        if(c.atoms()){
            if(structure_file==""){
                structure_file = s;
            } else {
                // Structure is present already
                // If this was set to topology file overwrite
                if(structure_file == top_file){
                    structure_file = s;
                } else if(!c.top()) {
                    throw Pteros_error("Only one structure file allowed!");
                }
            }

        }
    }

    if(traj_files.empty()) throw Pteros_error("At least one trajectory file is required!");

    // Print statistics of what we will read


    // Ensure we have tasks
    if(tasks.size()<1) throw Pteros_error("At least one task is required!");
    // Will read into the system of the first task
    System& system = tasks[0]->system;

    // To avoid reading top file twice
    if(structure_file==top_file) structure_file = "";

    if(structure_file!="" && top_file==""){
        // we have only structure but no topology
        system.load(structure_file);
    } else if(structure_file=="" && top_file!=""){
        // we have only topology but no structure
        system.load(top_file); // coordinates from top!
    } else if(structure_file!="" && top_file!=""){
        // we have both topology and structure
        system.load(structure_file);
        system.load(top_file); // No coordinates from top!
    } else {
        // No topology and no structure!
        // try using first TNG traj file as structure
        for(auto& s: traj_files){
            auto trj = Mol_file::recognize(s);
            auto c = trj->get_content_type();
            if(c.atoms() && c.traj()){
                structure_file = s;
                // We only need to load only atoms from TNG here
                log->debug("Using TNG trajectory file 's' to read structure...", s);
                trj->open('r');
                Frame fr;
                trj->read(&system, &fr, Mol_file_content().atoms(true));
                system.frame_append(fr);
                break;
            }
        }
        // If still no structure give up
        if(structure_file=="")
            Pteros_error("Structure AND/OR topology file is required!");
    }

    // Analysing which kind of tasks we have

    is_parallel = false;
    for(auto& task: tasks){
        if(task->is_parallel()){
            if(tasks.size()>1) throw Pteros_error("No other tasks can run if parallel task is present!");
            is_parallel = true;
            break;
        }
    }

    //-----------------------------------------
    // Actual processing starts here
    //-----------------------------------------    

    // Set buffer size
    int buf_size = options("buffer","10").as_int();    

    // Channel for frames
    Data_channel_ptr reader_channel(new Data_channel);
    reader_channel->set_buffer_size(buf_size);

    int Nproc = std::thread::hardware_concurrency();
    log->debug("Physical cores: {}", Nproc);
    log->debug("\tFile reading thread: 1");

    // Create traj file reader
    Traj_file_reader reader(options);
    // Start reader thread
    reader.run(traj_files, reader_channel);

    // Data container
    using Data_container_ptr = std::shared_ptr<Data_container>;
    Data_container_ptr data;

    // Processing depends on which tasks we have
    if(is_parallel){
        /* Single parallel task present
         * We run Nproc threads each with an instance of this task
         * and dispatch frames async to them from the reader channel
         * At the end we call Trajectory_reader::post_process()
         * to collect the results.
         * post_process() of individual instances are also called but
         * they only finalize particular instance.
         */

        // Start instances

        // We have Nproc-2 remote threads + this thread = Nproc-1 in total        
        int num_threads = Nproc-1;

        log->debug("\tThreads running parallel task: {}", num_threads+1);
        log->debug("\t({} separate + 1 master)", num_threads);

        // We have to reserve memory for all tasks in advance!
        // Otherwise due to reallocation of array pointers sent to threads may become invalid
        // which leads to f*cking misterious crashes!
        tasks.reserve(num_threads+1);

        // !!! Ensures that first valid frame is consumed by task[0] !!!
        // !!! Ensures that pre_process() is only called by task[0] on first valid frame !!!
        //
        // Jump remover (and probably some other things) have to initialize
        // itself on the first consumed frame
        // but each instance starts from unpredictable frame.
        // Thus we force task[0] to consume the first frame and call pre_process()
        // than we clone it with initialized state for other threads        
        tasks[0]->set_id(0);
        tasks[0]->driver->set_data_channel(reader_channel);
        tasks[0]->driver->init_with_first_frame(system);

        // Now spawn workers
        for(int i=1; i<=num_threads; ++i){ // task 0 will run in master thread, so start from 1
            // Clone provided task to make new independent instance
            tasks.push_back( Task_ptr(tasks[0]->clone()) );
            //cout << "clonned " << i << endl;
            tasks[i]->set_id(i);            
            tasks[i]->driver->set_data_channel(reader_channel);
            tasks[i]->driver->process_until_end_in_thread();            
        }

        // Process frame 0 which we got before spawning threads
        tasks[0]->driver->process_first_frame();

        // Run one worker in current thread        
        tasks[0]->driver->process_until_end();

        // Join all worker threads
        if(tasks.size()>1)
            for(int i=1; i<=num_threads; ++i) tasks[i]->driver->join_thread();


        // See which instance processed the last frame
        Frame_info last_info;
        int last_fr = -1;
        for(auto& t: tasks){
            if(t->driver->get_last_info().valid_frame > last_fr){
                last_fr = t->driver->get_last_info().valid_frame;
                last_info = t->driver->get_last_info();
            }
        }

        // Now collect results from all instances that consumed some frames
        vector<Task_ptr> resultive_tasks;
        for(int i=1; i<tasks.size();++i)
            if(tasks[i]->n_consumed) resultive_tasks.push_back(tasks[i]);

        log->debug("Collecting results from {} task instances...", resultive_tasks.size()+1);
        tasks[0]->collect_data(resultive_tasks);

        // Call post-process on master instance
        tasks[0]->post_process(last_info);


    } else {
        /* Only serial tasks are present
         * We make individual channels for each worker and feed the same frame
         * to each worker sequensially.
         * Each worker still runs in it's own thread.
         */        

        vector<Data_channel_ptr> worker_channels;

        if(tasks.size() > 1){
            // More than 1 consumer, start all of them in separate threads
            // Master thread will work as dispatcher

            log->debug("\tRunning {} serial tasks in separate threads", tasks.size());
            log->debug("\t(master thread is dispatching frames)");

            // We have to reserve memory for all channels in advance!
            // Otherwise due to reallocation of array pointers sent to threads may become invalid
            // which leads to f*cking misterious crashes!
            worker_channels.reserve(tasks.size());

            for(int i=0; i<tasks.size(); ++i){
                // Create new channel
                auto channel=std::make_shared<Data_channel>();
                channel->set_buffer_size(buf_size);
                worker_channels.push_back(channel);

                // Spawn worker
                tasks[i]->set_id(i);
                tasks[i]->driver->set_data_channel(worker_channels[i]);
                tasks[i]->driver->process_all_in_thread(system);
            }

            // Recieve all frames for reader channel and dispatch them to workers
            while(reader_channel->recieve(data)){
                for(auto &ch: worker_channels){
                    ch->send(data);
                }
            }

            // No more new frames, send stop to all workers
            for(auto &ch: worker_channels){
                ch->send_stop();
            }

            // Join all workers
            for(auto& t: tasks) t->driver->join_thread();

        } else {            
            // There is only one consumer, no need for multiple threads
            log->debug("\tRunning single serial task in master thread");
            tasks[0]->set_id(0);
            tasks[0]->driver->set_data_channel(reader_channel);
            tasks[0]->driver->process_all(system);
        }
    } // Dispatching frames

    // Join reader thread
    reader.join();

    log->debug("Trajectory processing finished!");

    auto end = chrono::steady_clock::now();

    log->info("Processing wall time: {}s", chrono::duration<double>(end-start).count() );

    // Print statistics
    if( is_parallel ){
        log->info("Number of frames processed by parallel task instances:");
        int tot = 0;
        for(int i=0; i<tasks.size(); ++i){
            log->info("\tInstance #{}: {}", i, tasks[i]->n_consumed);
            tot += tasks[i]->n_consumed;
        }
        log->info("\tTotal: {}", tot);
    } else {
        log->info("Number of frames processed by serial tasks:");
        for(int i=0; i<tasks.size(); ++i){
            log->info("\tTask #{}: {}", i,tasks[i]->n_consumed);
        }
    }
}


