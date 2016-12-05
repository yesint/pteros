/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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

#include <thread>
#include <boost/algorithm/string.hpp> // For to_lower
#include <boost/lexical_cast.hpp>

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


bool Trajectory_reader::is_frame_valid(int fr, float t){
    if(//------------------------------------------
        (first_frame<0 || fr>=first_frame)
        &&
        (first_time<0 || t>=first_time)
        &&
        (skip<0 || fr%skip==0)
      ){//------------------------------------------
        // This frame is valid
        return true;
    } else {
        // This frame is invalid
        return false;
    }
}

bool Trajectory_reader::is_end_of_interval(int fr, float t){
    if(//------------------------------------------
        (fr>last_frame && last_frame>=0)
            ||
        (t>last_time && last_time>=0)
      ){//------------------------------------------
        // End reached
        return true;
    } else {
        return false;
    }
}



Frame_info Trajectory_reader::dispatch_frames_to_task(const Task_ptr& task,
                                                const Data_channel_ptr& channel,
                                                const System& sys){

    std::shared_ptr<Data_container> data;

    bool pre_process_done = false;

    while(channel->recieve(data)){
        // put system only makes a copy of sys if it is not yet set in the task
        // So for the task[0] no excessive copy is made
        if(!pre_process_done) task->put_system(sys);
        task->put_frame(data->frame);
        if(!pre_process_done){
            task->pre_process();
            pre_process_done = true;
        }
        task->process_frame(data->frame_info);
    }

    // Call post_process
    task->post_process(data->frame_info);

    // Return last frame info for possible use in collector
    return data->frame_info;
}


void Trajectory_reader::process_value_with_suffix(const string& s, int* intval, float* floatval){
    size_t pos = s.find_last_of("0123456789");
    if(pos==string::npos) throw Pteros_error("A number with optional suffix required!");
    string val = s.substr(0,pos+1);
    string suffix = s.substr(pos+1);
    boost::algorithm::to_lower(val);
    boost::algorithm::to_lower(suffix);
    // Now analyze suffix
    if(intval!=nullptr && (suffix=="fr" || suffix=="")){
        *intval = boost::lexical_cast<int>(val);
        if(floatval) *floatval = -1.0;
    } else if(floatval!=nullptr) {
        if(intval) *intval = -1;
        *floatval = boost::lexical_cast<float>(val);
        if(suffix=="ps" || suffix=="t"){
            *floatval *= 1.0;
        } else if(suffix=="ns"){
            *floatval *= 1000.0;
        } else if(suffix=="us"){
            *floatval *= 1000000.0;
        } else if(suffix=="ms"){
            *floatval *= 1000000000.0;
        }
    } else if(floatval==nullptr && intval==nullptr){
        throw Pteros_error("Both int and float vals are NULL! WTF?");
    }
}

Trajectory_reader::Trajectory_reader(): collector(nullptr)
{

}

Trajectory_reader::Trajectory_reader(const Options &opt): options(opt), collector(nullptr)
{

}

void Trajectory_reader::run(){
    // Preparation stage

    cout << "Starting trajectory processing..." << endl;
    auto start = chrono::steady_clock::now();

    // Get files to work with
    auto file_list = options("f").as_strings();

    // Get all string parameters and find out how many different files we have
    // and arrange them structure->topology->traj
    string top_file = "";
    string structure_file = "";    

    for(string& s: file_list){                
        auto h = Mol_file::recognize(s);
        auto c = h->get_content_type();

        // traj file is always added as traj even if this is TNG
        if(c & MFC_TRAJ){
            traj_files.push_back(s);
            continue; // Avoid adding TNG twice also as structure file
        }

        if(c & MFC_TOP){
            if(top_file=="")
                top_file = s;
            else
                throw Pteros_error("Only one topology file allowed!");
        }

        if(c & MFC_ATOMS){
            if(structure_file==""){
                structure_file = s;
            } else {
                // Structure is present already
                // If this was set to topology file overwrite
                if(structure_file == top_file){
                    structure_file = s;
                } else if(!(c & MFC_TOP)) {
                    throw Pteros_error("Only one structure file allowed!");
                }
            }

        }
    }

    if(traj_files.empty()) throw Pteros_error("At least one trajectory file is required!");


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
            if(c & MFC_ATOMS && c & MFC_TRAJ){
                structure_file = s;
                // We only need to load only atoms from TNG here
                cout << "Using trajectory file " << s << " to read structure..." << endl;
                trj->open('r');
                Frame fr;
                trj->read(&system, &fr, MFC_ATOMS);
                system.frame_append(fr);
                break;
            }
        }
        // If still no structure give up
        if(structure_file=="")
            Pteros_error("Structure AND/OR topology file is required!");
    }


    // Get parameters for reading frames
    process_value_with_suffix(options("b","-1").as_string(),
                              &first_frame, &first_time);
    process_value_with_suffix(options("e","-1").as_string(),
                              &last_frame, &last_time);
    // Skip interval
    skip = options("skip","-1").as_int();

    // Check for custom start time and dt
    process_value_with_suffix(options("t0","-1").as_string(),
                              nullptr, &custom_start_time);
    process_value_with_suffix(options("dt","-1").as_string(),
                              nullptr, &custom_dt);
    if(custom_start_time>=0 && custom_dt==-1) custom_dt = 1;
    if(custom_start_time==-1 && custom_dt>=0) custom_start_time = 0;

    // Check if the range is valid
    if(first_frame>=0 && last_frame>=0 && last_frame<first_frame)
        throw Pteros_error("Last frame") << last_frame << " is smaller that first frame " << first_frame;
    if(first_time>=0 && last_time>=0 && last_time<first_time)
        throw Pteros_error("Last time") << last_time<< " is smaller that first time" << first_time;

    log_interval = options("log","-1").as_int();


    // Analysing which kind of tasks we have    


    is_parallel = false;
    for(auto& task: tasks){
        if(task->is_parallel()){
            if(tasks.size()>1) throw Pteros_error("No other tasks can run if parallel task is present!");
            is_parallel = true;
            break;
        }
    }

    if(is_parallel && !collector) throw Pteros_error("No collector function registered for parallel task!");


    //-----------------------------------------
    // Actual processing starts here
    //-----------------------------------------


    // Set buffer size
    int buf_size = options("buffer","10").as_int();    

    // Channel for frames
    Data_channel_ptr reader_channel(new Data_channel);
    reader_channel->set_buffer_size(buf_size);

    int Nproc = std::thread::hardware_concurrency();
    cout << "Physical cores: " << Nproc << endl;
    cout << "Threads used:" << endl;
    cout << "\tFile reading thread: 1" << endl;

    // Start reader thread    
    std::thread reader_thread( &Trajectory_reader::reader_thread_body, this, ref(reader_channel) );

    // Data container
    typedef std::shared_ptr<Data_container> Data_container_ptr;
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

        vector<std::thread> worker_threads;

        // Start instances

        // We have Nproc-2 remote threads + this thread = Nproc-1 in total        
        int num_threads = Nproc-1;

        cout << "\tThreads running parallel task: " << num_threads+1 << endl;
        cout << "\t(" << num_threads << " separate + 1 master)" << endl;

        // We have to reserve memory for all tasks in advance!
        // Otherwise due to reallocation of array pointers sent to threads may become invalid
        // which leads to f*cking misterious crashes!
        tasks.reserve(num_threads+1);

        for(int i=1; i<=num_threads; ++i){ // task 0 will run in master thread, so start from 1
            // Clone provided task to make new independent instance
            tasks.push_back( Task_ptr(tasks[0]->clone()) );

            worker_threads.push_back(
                        std::thread(
                                &Trajectory_reader::dispatch_frames_to_task,
                                this,
                                ref(tasks[i]),
                                ref(reader_channel),
                                ref(system)
                            )
                        );
        }


        // Run one worker in current thread
        Frame_info last_info = dispatch_frames_to_task(tasks[0],reader_channel,system);

        // Join all workers
        if(worker_threads.size()>0)
            for(auto& t: worker_threads) t.join();

        // Now collect results from all instances
        collector(last_info, tasks);


    } else {
        /* Only serial tasks are present
         * We make individual channels for each worker and feed the same frame
         * to each worker sequensially.
         * Each worker still runs in it's own thread.
         */        

        vector<std::thread> worker_threads;
        vector<Data_channel_ptr> worker_channels;

        if(tasks.size() > 1){
            // More than 1 consumer, start all of them in separate threads
            // Master thread will work as dispatcher

            cout << "\tRunning " << tasks.size() << " serial tasks in separate threads" << endl;
            cout << "\t(master thread is dispatcjing frames)" << endl;

            // We have to reserve memory for all channels in advance!
            // Otherwise due to reallocation of array pointers sent to threads may become invalid
            // which leads to f*cking misterious crashes!
            worker_channels.reserve(tasks.size());

            for(int i=0; i<tasks.size(); ++i){
                // Create new channel
                worker_channels.push_back(Data_channel_ptr(new Data_channel));
                // Set buffer size for this channel
                worker_channels.back()->set_buffer_size(buf_size);
                // Spawn thread
                worker_threads.push_back(
                            std::thread(
                                    &Trajectory_reader::dispatch_frames_to_task,
                                    this,
                                    ref(tasks[i]),
                                    ref(worker_channels[i]),
                                    ref(system)
                                )
                            );
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
            for(auto& t: worker_threads) t.join();

        } else {
            // There is only one consumer, no need for multiple threads
            cout << "\tRunning single serial task in master thread" << endl;
            dispatch_frames_to_task(tasks[0],reader_channel,system);
        }
    } // Dispatching frames

    // Join reader thread
    reader_thread.join();    

    cout << "Trajectory processing finished!" << endl;

    auto end = chrono::steady_clock::now();

    cout << "Processing time: " << chrono::duration<double>(end-start).count() << " s" << endl;
}

void Trajectory_reader::reader_thread_body(const Data_channel_ptr &channel){
    try {        
        int abs_frame = -1;
        int valid_frame = -1;
        // Saved first frame and time
        int saved_first_frame = -1;
        float saved_first_time = -1.0;

        bool finished = false;

        for(string& fname: traj_files){
            cout << "==> Reading trajectory " << fname << endl;

            auto trj = Mol_file::open(fname,'r');

            // Main loop over trajectory frames
            while(true){
                // To avoid excessive copy operations we allocate a shared pointer
                // and will load data into its storage
                std::shared_ptr<Data_container> data(new Data_container);

                // Load data to this container
                bool good = trj->read(NULL, &data->frame, MFC_TRAJ);

                // Check if EOF reached in trajectory
                if(!good) break;

                ++abs_frame; // Next absolute frame

                if(log_interval>0 && abs_frame%log_interval==0)
                    cout << "At frame " << abs_frame << endl;

                // If time stamps are overriden, use overrides
                if(custom_dt>=0){
                    data->frame.time = custom_start_time + custom_dt*abs_frame;
                }

                // Check if end of requested interval is reached
                if( is_end_of_interval(abs_frame,data->frame.time) ){
                    // Send stop to the queue
                    channel->send_stop();
                    finished = true;
                    // exit loop
                    break;
                }

                // Check if new frame falls into needed range of time.
                // If not go to next frame
                if( !is_frame_valid(abs_frame,data->frame.time) ) continue;

                // This is valid frame
                ++valid_frame;

                if(valid_frame==0){
                    // This is the very first valid frame, set start time
                    saved_first_frame = abs_frame;
                    saved_first_time = data->frame.time;
                }

                // Fill data container, which will be sent to the queue
                data->frame_info.absolute_time = data->frame.time;
                data->frame_info.absolute_frame = abs_frame;
                data->frame_info.valid_frame = valid_frame;                
                data->frame_info.first_frame = saved_first_frame;
                data->frame_info.first_time = saved_first_time;
                data->frame_info.last_frame = abs_frame;
                data->frame_info.last_time = data->frame.time;

                // Send frame to the queue
                channel->send(data);
            } // Over frames

            cout << "==> trajectory done" << endl;

            // If end reached break here too
            if(finished) break;

        } // Over trajectories

        // Send stop at the end
        channel->send_stop();

    } catch(Pteros_error e) {
        // Send stop if exception raised
        channel->send_stop();
        cout << "(ERROR) Execution of the reading thread stopped due to exception:" << endl;
        cout << e.what() << endl;
    }
}

