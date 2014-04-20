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
#include "pteros/analysis/trajectory_processor.h"
#include "pteros/core/pteros_error.h"
#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include "pteros/analysis/options_parser.h"
#include "pteros/core/format_recognition.h"
#include "pteros/core/mol_file.h"

using namespace pteros;
using namespace std;

Trajectory_processor::~Trajectory_processor(){
}

string Trajectory_processor::help(){
    return  "Trajectory processing options:\n"
            "General usage:\n"
            "\t--trajectory[ filename1 filename2 ... <processing options>]\n"
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
            "\t--first_frame <fr>\n"
            "\t\tfirst frame to read, default: 0\n"

            "\t--last_frame <fr>\n"
            "\t\tlast frame to read, default: -1 (up to the end)\n"

            "\t--first_time <t>\n"
            "\t\tfirst time step to read, ps, default: 0.0\n"

            "\t--last_time <t>\n"
            "\t\tlast time step to read, ps, default: -1 (up to the end)\n"

            "\t--window_size_frames sz:\n"
            "\t\tprocess by windows of size sz determined by frame.\n"

            "\t--window_size_time t:\n"
            "\t\tprocess by windows of size t determined by time.\n"

            "\t--skip <n>\n"
            "\t\tProcess only each n'th frame, default: -1 (process each frame)\n"

            "\t--custom_start_time <t>\n"
            "\t\tUse t as a starting, default: -1 (use value from first trajectory frame)\n"
            "\t\tUseful if trajectory does not contain time stamps\n"
            "\t\tor if the starting time is incorrect.\n"
            "\t\tIf set and custom_dt is not given sets custom_dt to 1.0!\n"

            "\t--custom_dt <t>\n"
            "\t\tUsed as a time step, default: -1 (use value from trajectory)\n"
            "\t\tUseful if trajectory does not contain time stamps.\n"
            "\t\tIf set and custom_start_time is not given sets custom_start_time to 0.0!\n"

            "\t--log_interval n\n"
            "\t\tPrints logging information each n frames, default -1 (no logging)\n"

            "\t--dump_input filename\n"
            "\t\tDumps input in JSON format to specified file\n";
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

void Trajectory_processor::run(){
    cout << "Starting trajectory processing..." << endl;

    // See if thre are some connected consumers
    if(consumers.size()==0) throw Pteros_error("No consumers are connected to trajectory processor!");
    cout << "Connected " << consumers.size() << " consumers" << endl;

    // If asked to dump input in json format, do it
    string dump_file = options->get_value<string>("dump_input","");
    if(dump_file!=""){
        cout << "Dumping input to " << dump_file << "..." << endl;
        ofstream f(dump_file.c_str());
        f << options->to_json_string() << endl;
        f.close();
    }

    // Work with trajectory block
    // The block must include structre file,
    // one or more trajectory files and optional topology file
    // Get trajectory node
    Options_tree* trj = &options->get_option("trajectory");

    // Get all string parameters and find out how many different files we have
    // and arrange them structure->topology->traj
    string top_file = "";
    string structure_file = "";
    traj_files.clear();
    for(string s: trj->get_values<string>("")){
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
    window_size_frames = trj->get_value<int>("window_size_frames",-1);
    window_size_time = trj->get_value<float>("window_size_time",-1);
    // Determine range for this group
    first_frame = trj->get_value<int>("first_frame",-1);
    last_frame = trj->get_value<int>("last_frame",-1);
    first_time = trj->get_value<float>("first_time",-1);
    last_time = trj->get_value<float>("last_time",-1);
    skip = trj->get_value<int>("skip",-1);
    // Check for custom start time and dt
    custom_start_time = trj->get_value<float>("custom_start_time",-1);
    custom_dt = trj->get_value<float>("custom_dt",-1);
    if(custom_start_time>=0 && custom_dt==-1) custom_dt = 1;
    if(custom_start_time==-1 && custom_dt>=0) custom_start_time = 0;

    // Check if the range is valid
    if(first_frame>=0 && last_frame>=0 && last_frame<first_frame)
        throw Pteros_error("Last frame") << last_frame << " is smaller that first frame " << first_frame;
    if(first_time>=0 && last_time>=0 && last_time<first_time)
        throw Pteros_error("Last time") << last_time<< " is smaller that first time" << first_time;

    log_interval = trj->get_value<int>("log_interval",-1);

    //-----------------------------------------
    // Actual processing starts here
    //-----------------------------------------
    typedef std::shared_ptr<Data_channel> Data_channel_ptr;

    // Set buffer size
    int buf_size = options->get_value<int>("buffer_size",10);
    cout << "Using frame buffers of size " << buf_size << endl;
    channel.set_buffer_size(buf_size);

    // Start reader thread
    boost::thread reader_thread( boost::bind(&Trajectory_processor::reader_thread_body,this) );

    boost::thread_group worker_threads;
    vector<Data_channel_ptr> worker_channels;

    if(consumers.size() > 1){
        // More than 1 consumer, start all of them in separate threads
        for(int i=0; i<consumers.size(); ++i){
            // Create channel
            worker_channels.push_back(Data_channel_ptr(new Data_channel));
            // Set buffer size for this consumer
            worker_channels.back()->set_buffer_size(buf_size);
            // Spawn thread
            worker_threads.create_thread(
                        boost::bind(&Consumer_base::run_in_thread,
                                              consumers[i],
                                              worker_channels[i]
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

        //!! Important !!
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
        worker_threads.join_all();
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
