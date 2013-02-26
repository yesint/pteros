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

#include <fstream>
#include "pteros/analysis/trajectory_processor.h"
#include "pteros/core/pteros_error.h"
#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include "pteros/analysis/options_parser.h"
#include "pteros/core/format_recognition.h"
#include "pteros/core/mol_file.h"

using namespace pteros;
using namespace std;

void Trajectory_processor::init(){

}

Trajectory_processor::~Trajectory_processor(){
}

void Trajectory_processor::print_help(){
    cout << "Note about nested options:\n"
            "--------------------------\n"
            "Nested options could be used by putting arguments of the parent options\n"
            "into square brackets like this:\n"
            "--parent [ arg1 arg2 --nested1 nested_arg1 nested_arg2 ]"

            "General options:\n"
            "----------------\n"
            "--help:\n\tPrint this help message\n"

            "--json filename:\n\tRead options from specified JSON file.\n"
            "\n"

            "Options for trajectory processing:\n"
            "----------------------------------\n"
            "--trajectory filename1 filename2 ... <sub-options>:\n"
            "\tRequired. Group of files, which include:\n"
            "\t* structure file (PDB or GRO, required),\n"
            "\t* topology file (preprocessed Gromacs TOP, optional)\n"
            "\t* one or more trajectory files (TRR or XTC, required).\n"
            "\tFiles may appear in any order, but trajectory files will be processed\n"
            "\tin the order of their appearance.\n"
            "\tThis option could be given several times.\n"
            "\tThe following sub-options may appear inside --trajectory:\n\n"

            "\t--range [frame_range|time_range] begin end:\n"
            "\t\tOptional. Part of trajectory to process expressed in either frames or time in ps.\n"

            "\t--window [frame_window|time_window] sz:\n"
            "\t\tprocess by windows of size sz determined by frame of by time in ps.\n"

            "\t--name str:\n"
            "\t\tOptional. Used in the names of output files.\n"

            "\n"

            "--async [true|false]\n"
            "\tOptional. Asynchronous reading of trajectory frames. True by default.\n"

            "--log_interval n\n"
            "\tOptional. Print logging info each n frames.\n"

            "--dump_input file\n"
            "\tOptional. Dumps input in JSON format to specified file.\n"

            "\n"
         << endl;
}

bool Trajectory_processor::check_time_range(int fr, float t){
    if(//------------------------------------------
        (first_frame<0 || fr>=first_frame)
        &&
        (last_frame<0 || fr<=last_frame)
        &&
        (first_time<0 || t>=first_time)
        &&
        (last_time<0 || t<=last_time)
        &&
        (skip<0 || fr%skip==0)
      ){//------------------------------------------
        // This frame is valid
        return true;
    } else {
        // This frame is invalid. See if we go past the end
        if(//------------------------------------------
            (fr>last_frame && last_frame>=0)
                ||
            (t>last_time && last_time>=0)
          ){//------------------------------------------
            // End reached, set stop flag
            stop_requested = true;
        }
        // And rturn false in any case
        return false;
    }
}

void Trajectory_processor::read_single_trajectory(string& fname){
    bool term = false; // Termination flag
    try {
        cout << "==> Reading trajectory " << fname << endl;

        boost::shared_ptr<Mol_file> trj = io_factory(fname,'r');
        if(!trj->get_content_type().trajectory)
            throw Pteros_error("File "+fname+" is not a trajectory!");

        while(true){
            // To avoid eccessive copy operations we allocate a shared pointer
            // and will load data into its storage
            boost::shared_ptr<Data_container> data(new Data_container);

            // Load data to this container
            // This may safely happen outside the lock
            Mol_file_content c;
            c.trajectory = true;
            bool good = trj->read(NULL,&data->frame,c);            

            // Check if EOF reached in trajectory
            if(!good) break;

            {
                boost::mutex::scoped_lock lk(buffer_mutex);

                ++abs_frame; // Next absolute frame

                // If time stamps are overriden, use overrides
                if(custom_dt>=0){
                    data->frame.t = custom_start_time + custom_dt*abs_frame;
                }

                // Check if new frame falls into needed range of time.
                bool ok = check_time_range(abs_frame,data->frame.t);
                // check_time_range may set stop flag if the end of range is reached, so check
                if(stop_requested) {                    
                    break;
                }
                // If we are here, then stop is not requested, but frame may still be out of range
                if(!ok) continue; // Go to next frame

                // Wait until buffer size becomes smaller then buffer_size or until stop requested
                while(buffer.size()>=max_buffer_size && !stop_requested){
                    buffer_cond.wait(lk);
                }                
            }

            // If stop requested break here
            if(stop_requested){
                break;
            }

            // This is valid frame
            ++valid_frame;

            if(valid_frame==0){
                // This is the very first valid frame, set start time
                saved_first_frame = abs_frame;
                saved_first_time = data->frame.t;
            }                       

            // Fill data container, which will be sent to the queue                        
            data->frame_info.absolute_time = data->frame.t;
            data->frame_info.absolute_frame = abs_frame;
            data->frame_info.valid_frame = valid_frame;
            data->frame_info.win_size_frames = window_size_frames;
            data->frame_info.win_size_time = window_size_time;
            data->frame_info.first_frame = saved_first_frame;
            data->frame_info.first_time = saved_first_time;            

            // Block mutex and add new frame to the queue            
            {
                boost::mutex::scoped_lock lock(buffer_mutex);

                if(log_interval>0 && abs_frame%log_interval==0)
                    cout << "Loaded frame " << abs_frame << endl;

                buffer[valid_frame] = data;
                // Set access count for this frame to empty set (nobody asked for it yet)
                frame_access_count[valid_frame] = vector<bool>(consumers.size(),false);               

                // Notify that new frame arrived
                buffer_cond.notify_all();                
            }

        }

        // On exit from the reading loop set stop condition and notify any waiting consumers
        {
            boost::mutex::scoped_lock lock(buffer_mutex);
            stop_requested = true;
            buffer_cond.notify_all();
        }

        cout << "==> reading done" << endl;

    } catch(Pteros_error e) {
        cout << e.what() << endl;
    }
}

void Trajectory_processor::start_threads(vector<string>& fnames){

    cout << "--> There are " << consumers.size() << " consumers" << endl;

    stop_requested = false;
    max_buffer_size = options->get_value<int>("buffer_size",10);

    alive_consumers.resize(consumers.size());
    for(int i=0; i<consumers.size(); ++i) alive_consumers[i] = true;

    // Create consumer threads
    for(int i=0; i<consumers.size(); ++i){
        consumer_threads.create_thread( boost::bind(&Consumer::run,consumers[i]) );
    }

    // Read supplied trajectories one by one
    valid_frame = -1;
    abs_frame = -1;
    for(int i=0; i<fnames.size(); ++i){
        read_single_trajectory(fnames[i]);
    }

    // Wait for all consumer threads to finish
    consumer_threads.join_all();
}

// Provides frame #fr for the consumer #id
boost::shared_ptr<Data_container> Trajectory_processor::frame_provider(int fr, int id){
    // Lock mutex
    boost::mutex::scoped_lock lk(buffer_mutex);           

    // Wait until frame fr appears in the buffer or until stop requested
    while(buffer.count(fr)==0 && !stop_requested){
        buffer_cond.wait(lk); // This condition will be notified from the reading thread
    }

    // If stop requested and there is no such frame in the buffer, exit
    if(stop_requested && buffer.count(fr)==0){
        // We will not read from the buffer, so we need to allocate new return object
        boost::shared_ptr<Data_container> ret(new Data_container);
        // We fill window info and first time (this is frame-independent).
        ret->frame_info.win_size_frames = window_size_frames;
        ret->frame_info.win_size_time = window_size_time;
        ret->frame_info.first_frame = saved_first_frame;
        ret->frame_info.first_time = saved_first_time;
        ret->stop = true; // Signal for consumer to stop
        // there is still no absolute_time, because we have no frame to extract it from
        // this will be filled in the consumer itself
        return ret;
    }

    // Get frame from the queue and modify counter
    boost::shared_ptr<Data_container> ret =  buffer[fr];
    frame_access_count[fr][id] = true;

    // If all alive consumers accessed this frame, delete it from the buffer
    int a = 0; // Number of alive consumers
    int b = 0; // Number of accesses
    for(int i=0; i<consumers.size(); ++i){
        if(alive_consumers[i]){
            ++a;
            if(frame_access_count[fr][i]) ++b;
        }
    }
    if(a==b){        
        buffer.erase(fr);
        frame_access_count.erase(fr);
        // And notify if reader is waiting for empty slot in buffer
        buffer_cond.notify_all();
    }

    return ret; // send data to consumer
}

void Trajectory_processor::consumer_finished(int id){
    {
        boost::mutex::scoped_lock lk(buffer_mutex);
        alive_consumers[id] = false;

        cout << "Consumer " << id << " finished" << endl;

        // If any other consumer is still alive just exit
        int a = 0;
        for(int i=0; i<consumers.size(); ++i){
            if(alive_consumers[i]) ++a;
        }
        if(a>0) return;

        // If we are here, then all consumers are finished
        // Signal to reading thread to stop
        cout << "All consumers finished!" << endl;
        cout << "Remaining buffer size: " << buffer.size() << endl;
        // This will signal for reader thread to stop if it is still running
        stop_requested = true;
        // Notify reader to release it from waiting
        buffer_cond.notify_one();
    }
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

    log_interval = options->get_value<int>("log_interval",0);

    // Work with trajectory block
    // The block reads one trajectory, which must include structre file,
    // one or more trajectory files and optional topology file   
    // Get trajectory node
    Options_tree* trj = &options->get_option("trajectory");

    // Get all string parameters and find out how many different files we have
    // and arrange them structure->topology->traj
    string top_file = ""; // Only one allowed
    string structure_file = "";
    vector<string> traj_files;
    for(string s: trj->get_values<string>("")){
        switch(recognize_format(s)){
        case PDB_FILE:
        case GRO_FILE:
            if(structure_file!="") throw Pteros_error("Only one structure file is allowed!");
            structure_file = s;
            break;
        case TOP_FILE:
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

    // Now read everything in sequence

    // Read structure file
    if(structure_file=="") throw Pteros_error("Structure file is required!");
    // We will read this file into the system of first supplied consumer
    // and then copy that system to all other consumers
    System* sys1 = consumers[0]->get_system();
    sys1->clear();
    sys1->load(structure_file);

    cout << "Copying system data to consumers..." << endl;
    for(int i=1; i<consumers.size(); ++i)
        *(consumers[i]->get_system()) = *sys1;

    // Read topology if provided
    if(top_file!=""){
        simulation = boost::shared_ptr<Simulation>(new Simulation);
        // Load data into simulation using system from first consumer
        simulation->create(*(consumers[0]->get_system()),top_file);
        // Copy pointer to all consumers
        for(int i=0; i<consumers.size(); ++i){
            consumers[i]->set_simulation(simulation);
        }
    }

    // Now read trajectories
    if(traj_files.empty()) throw Pteros_error("At least one trajectory file is required!");

    // See if window processing is requested in this group
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

    // Run reading and processing threads
    start_threads(traj_files);

    cout << "Trajectory processing finished!" << endl;
}

void Trajectory_processor::add_consumer(Consumer *p){
    consumers.push_back(p);
    consumers.back()->set_id(consumers.size()-1);
}
