/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/





#include "pteros/analysis/trajectory_reader.h"
#include "message_channel.h"
#include "data_container.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/file_handler.h"
#include "task_driver.h"
#include "traj_file_reader.h"
#include "pteros/core/logging.h"
#include <thread>

using namespace pteros;
using namespace std;


string TrajectoryReader::help(){
    return
R"(Trajectory processing options:
General usage:
    -f filename1 filename2 ... <processing options>
Files:
    * Exactly one structure file (PDB or GRO)
      If not specified, topology TPR file or TNG trajectory file
      must be given instead.
    * Gromacs topology file (TPR)
      If structure file is also present only topology is read from this file.
      If structure file is not present the coordinates are also read.
    * One or more trajectory files (TRR, XTC, TNG or DCD).
      TNG files also contain the structure, so if no structure file
      is given the structure is read from the first TNG file.

    Files may appear in any order, but trajectory files will be processed
    in the order of their appearance.

Processing options:

    -path <string>
        optional path which will be prepended to all data files, default: empty string
    -b <value[suffix]>
        beginning of processing (starting frame or time), default: 0
    -e <value[suffix]>
        end of processing (end frame or time), default: -1 (up to the end)
    -skip <n>
        Process only each n'th frame, default: -1 (process each frame)
    -t0 <t[suffix]>
        Custom starting time, default: -1 (use value from first frame)
        Useful if trajectory does not contain time stamps
        or if the starting time is incorrect.
        If set and dt is not given sets dt to 1.0!
    -dt <t[suffix]>
        Cutom time step, default: -1 (use value from trajectory)
        Useful if trajectory does not contain time stamps.
        If set and start is not given sets start to 0.0!
    -log <n>
        Prints logging information on each n-th frame, default: -1 (no logging)
    -buffer <n>
        Number of frames, which are kept in memory, default: 10
        Only touch this if individual frames are very large.

Suffixes:
    All parameters marked as <value[suffix]> accept the following optional suffixes:
        (no suffix) - value is in frames
        fr - value is in frames
        t - value is time in picoseconds (value used as is)
        ps - value is time in picoseconds (value used as is)
        ns - value is time in nanoseconds (value multiplied by 10^3)
        us - value is time in microseconds (value multiplied by 10^6)
        ms - value is time in milliseconds (value multiplied by 10^9)
    Parameters marked as <t[suffix]> does not accept fr suffix.
    In this case no suffix means ps.
)";
}

void TrajectoryReader::add_task(TaskBase *task){
    tasks.emplace_back(task);
}

void TrajectoryReader::add_task(const Task_ptr &task){
    tasks.push_back(task);
}



TrajectoryReader::TrajectoryReader()
{

}

TrajectoryReader::TrajectoryReader(const Options &opt): options(opt)
{

}

void TrajectoryReader::set_options(const Options &opt){
    options = opt;
}

void TrajectoryReader::run(){    
    // Separate logger (not registered since only used here)
    auto log = create_logger("trj_reader");

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
        auto h = FileHandler::recognize(s);
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
                throw PterosError("Only one topology file allowed!");
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
                    throw PterosError("Only one structure file allowed!");
                }
            }

        }
    }

    if(traj_files.empty()) throw PterosError("At least one trajectory file is required!");

    // Ensure we have tasks
    if(tasks.size()<1) throw PterosError("At least one task is required!");
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
            auto trj = FileHandler::recognize(s);
            auto c = trj->get_content_type();
            if(c.atoms() && c.traj()){                
                structure_file = s;
                // We only need to load only atoms from TNG here
                log->debug("Using trajectory file '{}' to read structure...", s);
                trj->open('r');
                Frame fr;
                trj->read(&system, &fr, FileContent().atoms(true));
                system.frame_append(fr);
                break;
            }
        }

        // If still no structure give up
        if(structure_file=="") throw PterosError("Structure AND/OR topology file is required!");
    }

    // Analysing which kind of tasks we have

    is_parallel = false;
    for(auto& task: tasks){
        if(task->is_parallel()){
            if(tasks.size()>1) throw PterosError("No other tasks can run if parallel task is present!");
            is_parallel = true;
            break;
        }
    }

    // Print summary of files we are going to process
    if(log->level() <= spdlog::level::debug){
        log->debug("Summary of files to be processed:");
        log->debug("\tStructure file:\t'{}'",structure_file);
        log->debug("\tTopology file:\t'{}'",top_file);
        log->debug("\tTrajectory files:");
        for(auto& f: traj_files) log->debug("\t\t{}",f);
    }

    //-----------------------------------------
    // Actual processing starts here
    //-----------------------------------------    

    // Set buffer size
    int buf_size = options("buffer","10").as_int();    

    // Channel for frames
    DataChannel_ptr reader_channel(new DataChannel);
    reader_channel->set_buffer_size(buf_size);

    int Nproc = std::thread::hardware_concurrency();
    log->debug("Physical cores: {}", Nproc);
    log->debug("\tFile reading thread: 1");

    // Create traj file reader
    Traj_file_reader reader(options, system.num_atoms());
    // Start reader thread
    reader.run(traj_files, reader_channel);

    // Data container
    using Data_container_ptr = std::shared_ptr<DataContainer>;
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

        tasks[0]->set_id(0);
        tasks[0]->driver->set_data_channel_and_system(reader_channel,system);

        // Call user-defined init before spawning tasks. System is already set.
        // For parallel tasks jump remover is initialized inside this call
        // and then is cloned around
        tasks[0]->before_spawn_handler();

        // Now spawn other workers
        for(int i=1; i<=num_threads; ++i){ // task 0 will run in master thread, so start from 1
            // Clone provided task to make new independent instance
            tasks.emplace_back(tasks[0]->clone());
            tasks[i]->set_id(i);
            tasks[i]->driver->set_data_channel_and_system(reader_channel,system);
            tasks[i]->driver->process_until_end_in_thread();
        }

        // Run one worker in current thread        
        tasks[0]->driver->process_until_end();

        // Join all worker threads
        if(tasks.size()>1)
            for(int i=1; i<=num_threads; ++i) tasks[i]->driver->join_thread();


        // Now collect results from all instances that consumed some frames
        vector<Task_ptr> resultive_tasks;
        int n_total = 0;
        for(int i=1; i<tasks.size();++i){
            if(tasks[i]->n_consumed) resultive_tasks.push_back(tasks[i]);
            n_total += tasks[i]->n_consumed;
        }

        log->debug("Collecting results from {} task instances...", resultive_tasks.size()+1);
        tasks[0]->collect_data(resultive_tasks,n_total);

    } else {
        /* Only serial tasks are present
         * We make individual channels for each worker and feed the same frame
         * to each worker sequensially.
         * Each worker still runs in it's own thread.
         */        

        vector<DataChannel_ptr> worker_channels;

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
                auto channel=std::make_shared<DataChannel>();
                channel->set_buffer_size(buf_size);
                worker_channels.push_back(channel);

                // Configure worker
                tasks[i]->set_id(i);
                tasks[i]->driver->set_data_channel_and_system(worker_channels[i],system);
                tasks[i]->driver->process_until_end_in_thread();
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
            tasks[0]->driver->set_data_channel_and_system(reader_channel,system);
            tasks[0]->driver->process_until_end();
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






