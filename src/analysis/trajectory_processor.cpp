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

#include <fstream>
#include <thread>
#include <functional>

#include "pteros/analysis/trajectory_processor.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/format_recognition.h"
#include "pteros/core/mol_file.h"

#include <boost/algorithm/string.hpp> // For to_lower
#include <boost/lexical_cast.hpp>

using namespace pteros;
using namespace std;

Trajectory_processor::~Trajectory_processor(){
}

string Trajectory_processor::help(){
    return  "Trajectory processing options:\n"
            "General usage:\n"
            "\t-t filename1 filename2 ... <processing options>\n"
            "Files:\n"
            "\t* Exactly one structure file (PDB or GRO)\n"
            "\t  If not specified, topology PTTOP file must be given instead.\n"
            "\t* Topology PTTOP file (converted from Gromacs .tpr by tpr2pteros.py)\n"
            "\t  If structure file is also present only topology is read from this file.\n"
            "\t  If structure file is not present the coordinates are also read.\n"
            "\t* One or more trajectory files (TRR, XTC of DCD).\n"
            "\n"
            "\tFiles may appear in any order, but trajectory files will be processed\n"
            "\tin the order of their appearance.\n\n"

            "Processing options:\n"
            "\t-b <value[suffix]>\n"
            "\t\tbeginning of processing (starting frame or time), default: 0\n"

            "\t-e <value[suffix]>\n"
            "\t\tend of processing (end frame or time), default: -1 (up to the end)\n"

            "\t-w <value[suffix]>\n"
            "\t\tprocess by windows of given size in frames or time.\n"

            "\t-skip <n>\n"
            "\t\tProcess only each n'th frame, default: -1 (process each frame)\n"

            "\t-t0 <t[suffix]>\n"
            "\t\tCustom starting time, default: -1 (use value from first frame)\n"
            "\t\tUseful if trajectory does not contain time stamps\n"
            "\t\tor if the starting time is incorrect.\n"
            "\t\tIf set and dt is not given sets dt to 1.0!\n"

            "\t-dt <t[suffix]>\n"
            "\t\tCutom time stap, default: -1 (use value from trajectory)\n"
            "\t\tUseful if trajectory does not contain time stamps.\n"
            "\t\tIf set and start is not given sets start to 0.0!\n"

            "\t-log <n>\n"
            "\t\tPrints logging information on each n-th frame, default: -1 (no logging)\n"

            "\t-buffer <n>\n"
            "\t\tNumber of frames, which are kept in memory, default: 10\n"
            "\t\tOnly touch this if individual frames are very large.\n"

            "Suffixes:\n"
            "\tAll parameters markes as <value[suffix]> accept following optional suffixes:\n"
            "\t\t(no suffix) - value is in frames\n"
            "\t\tfr - value is in frames\n"
            "\t\tt - value is time in picoseconds (value used as is)\n"
            "\t\tps - value is time in picoseconds (value used as is)\n"
            "\t\tns - value is time in nanoseconds (value multiplied by 10^3)\n"
            "\t\tus - value is time in microseconds (value multiplied by 10^6)\n"
            "\t\tms - value is time in milliseconds (value multiplied by 10^9)\n"
            "\tParameters markes as <t[suffix]> does not accept fr suffix.\n"
            "\tIn this case no suffix means ps.\n"
            ;
}

bool Trajectory_processor::is_frame_valid(int fr, float t){
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

bool Trajectory_processor::is_end_of_interval(int fr, float t){
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

void Trajectory_processor::add_consumer(Consumer_base *p){
    consumers.push_back(p);
    consumers.back()->set_id(consumers.size()-1);
}

void process_value_with_suffix(const string& s, int* intval, float* floatval){
    int pos = s.find_last_of("0123456789");
    if(pos==string::npos) throw Pteros_error("A number with optional suffix required!");
    string val = s.substr(0,pos+1);
    string suffix = s.substr(pos+1);
    boost::algorithm::to_lower(val);
    boost::algorithm::to_lower(suffix);
    // Now analyze suffix
    if(intval!=nullptr && (suffix=="fr" || suffix=="")){
        *intval = boost::lexical_cast<int>(val);
        *floatval = -1.0;
    } else if(floatval!=nullptr) {
        *intval = -1;
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

void Trajectory_processor::run(){
    cout << "Starting trajectory processing..." << endl;

    // See if thre are some connected consumers
    if(consumers.size()==0) throw Pteros_error("No consumers are connected to trajectory processor!");
    cout << "Connected " << consumers.size() << " consumers" << endl;    

    // Get files to work with
    auto file_list = options("f").as_strings();

    // Get all string parameters and find out how many different files we have
    // and arrange them structure->topology->traj
    string top_file = "";
    string structure_file = "";    
    for(string& s: file_list){
        switch(recognize_format(s)){
        case PDB_FILE:
        case GRO_FILE:
            if(structure_file!="") throw Pteros_error("Only one structure file is allowed!");
            structure_file = s;
            break;
        case PTTOP_FILE:
            if(top_file!="") throw Pteros_error("Only one topology file is allowed!");
            top_file = s;
            break;
        case TRR_FILE:
        case XTC_FILE:
        case DCD_FILE:
            traj_files.push_back(s);
            break;
        }
    }

    if(traj_files.empty()) throw Pteros_error("At least one trajectory file is required!");

    // We will read into the system of first supplied consumer
    // and then copy that system to all other consumers
    System* sys1 = consumers[0]->get_system();
    sys1->clear();

    if(structure_file!="" && top_file==""){
        // we have only structure but no topology
        sys1->load(structure_file);
    } else if(structure_file=="" && top_file!=""){
        // we have only topology but no structure
        sys1->load(top_file); // coordinates and stucture from top!
    } else if(structure_file!="" && top_file!=""){
        // we have both topology and structure
        sys1->load(structure_file);
        sys1->load(top_file); // No coordinates and stucture from top!
    } else {
        Pteros_error("Structure AND/OR topology file is required!");
    }

    // Copy system to other consumers if needed
    if(consumers.size()>1){
        cout << "Copying system data to consumers..." << endl;
        for(int i=1; i<consumers.size(); ++i){
            *(consumers[i]->get_system()) = *sys1; // deep copying
        }
    }

    // Get parameters
    // See if window processing is requested    
    process_value_with_suffix(options("w","-1").as_string(),
                              &window_size_frames, &window_size_time);

    // Determine range for this group
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

    //-----------------------------------------
    // Actual processing starts here
    //-----------------------------------------
    typedef std::shared_ptr<Data_channel> Data_channel_ptr;

    // Set buffer size
    int buf_size = options("buffer","10").as_int();
    cout << "Using frame buffers of size " << buf_size << endl;

    channel.set_buffer_size(buf_size);

    // Start reader thread    
    std::thread reader_thread( &Trajectory_processor::reader_thread_body, this );

    vector<std::thread> worker_threads;
    vector<Data_channel_ptr> worker_channels;

    if(consumers.size() > 1){
        // More than 1 consumer, start all of them in separate threads
        for(int i=0; i<consumers.size(); ++i){
            // Create channel
            worker_channels.push_back(Data_channel_ptr(new Data_channel));
            // Set buffer size for this consumer
            worker_channels.back()->set_buffer_size(buf_size);
            // Spawn thread            
            worker_threads.push_back(
                        std::thread(
                            // We need bind here - doesn't work without it
                            std::bind(&Consumer_base::run_in_thread,
                                              consumers[i],
                                              worker_channels[i]
                                              )
                            )
                        );
        }

        // Now recieve frames from the queue until reader sends stop
        std::shared_ptr<Data_container> data;
        while(channel.recieve(data)){
            for(auto &ch: worker_channels){
                ch->send(data);
            }
        }

        // If we are here than reader thread sent a stop to the queue
        // Consume all remaining frame
        while(!channel.empty()){
            channel.recieve(data);
            for(auto &ch: worker_channels){
                ch->send(data);
            }
        }

        // No more new frames, send stop to all consumers
        for(auto &ch: worker_channels){
            ch->send_stop();
        }
    } else {
        // There is only one consumer, no need for multiple threads
        // Run pre-process

        // Important !!
        // Try block here does not catch exceptions inside pre_process
        // because this could be called from externally loaded compiled plugin!
        // errors should be catched in the Consumer itself!
        consumers[0]->pre_process_handler();

        std::shared_ptr<Data_container> data;
        while(channel.recieve(data)){
            consumers[0]->consume_frame(data);
        }
        while(!channel.empty()){
            channel.recieve(data);
            consumers[0]->consume_frame(data);
        }
        // Run post-process with last supplied data
        consumers[0]->post_process_handler(data->frame_info);
    }

    // Join all threads
    reader_thread.join();    
    if(consumers.size() > 1){        
        for(auto& t: worker_threads) t.join();
    }

    cout << "Trajectory processing finished!" << endl;
}

void Trajectory_processor::reader_thread_body(){
    try {
        // Reading content flag
        Mol_file_content content;
        content.trajectory = true;

        int abs_frame = -1;
        int valid_frame = -1;
        // Saved first frame and time
        int saved_first_frame = -1;
        float saved_first_time = -1.0;

        bool finished = false;

        for(string& fname: traj_files){
            cout << "==> Reading trajectory " << fname << endl;

            auto trj = io_factory(fname,'r');

            // Main loop over trajectory frames
            while(true){
                // To avoid eccessive copy operations we allocate a shared pointer
                // and will load data into its storage
                std::shared_ptr<Data_container> data(new Data_container);

                // Load data to this container
                bool good = trj->read(NULL,&data->frame,content);

                // Check if EOF reached in trajectory
                if(!good) break;

                ++abs_frame; // Next absolute frame

                if(log_interval>0 && abs_frame%log_interval==0)
                    cout << "Loaded frame " << abs_frame << endl;

                // If time stamps are overriden, use overrides
                if(custom_dt>=0){
                    data->frame.time = custom_start_time + custom_dt*abs_frame;
                }

                // Check if end of requested interval is reached
                if( is_end_of_interval(abs_frame,data->frame.time) ){
                    // Send stop to the queue
                    channel.send_stop();
                    finished = true;
                    // End exit loop
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
                data->frame_info.win_size_frames = window_size_frames;
                data->frame_info.win_size_time = window_size_time;
                data->frame_info.first_frame = saved_first_frame;
                data->frame_info.first_time = saved_first_time;
                data->frame_info.last_frame = abs_frame;
                data->frame_info.last_time = data->frame.time;

                // There is no data about window starts and ends here!
                // This information is managed by the Consumer!

                // Send frame to the queue
                channel.send(data);                
            } // Over frames

            cout << "==> reading done" << endl;

            // If end reached break here too
            if(finished) break;

        } // Over trajectories

        // Send stop at the end
        channel.send_stop();

    } catch(Pteros_error e) {
        // Send stop if exception raised
        channel.send_stop();        
        cout << "(ERROR) Execution of the reading thread stopped due to exception!" << endl;
        cout << e.what() << endl;
    }
}
