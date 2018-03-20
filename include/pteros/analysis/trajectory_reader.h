/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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



#ifndef TRAJECTORY_READER_H
#define TRAJECTORY_READER_H

#include <string>
#include <functional>
#include "pteros/analysis/options.h"
#include "pteros/analysis/task_base.h"


namespace pteros {

using Task_ptr = std::shared_ptr<Task_base> ;
using Data_channel = Message_channel<std::shared_ptr<pteros::Data_container> > ;
using Data_channel_ptr = std::shared_ptr<Data_channel> ;

/** The base class for trajectory processing
*   It provides facilities for loading large trajectories by frames
*   and to analyze each frame by user-defined function.
*   The range of processing could be given
*   by frame number or by physical time.
*/
class Trajectory_reader {
public:

        /// Default constructor
        Trajectory_reader();
        /// Constructor with options
        Trajectory_reader(const Options& opt);
        /// Destructor
        virtual ~Trajectory_reader(){}

        /// Pass options
        void set_options(const Options& opt){
            options = opt;
        }

        /// Read trajectory
        virtual void run();

        /// Print summary of allowed options
        std::string help();

        /// Adds new task
        void add_task(Task_base* task){
            tasks.emplace_back(task);
        }

        void add_task(const Task_ptr& task){
            tasks.push_back( task );
        }        

private:

        // Options
        Options options;

        //void reader_thread_body(const Data_channel_ptr &channel);

        std::vector<std::string> traj_files;        

        std::vector<Task_ptr> tasks;

        bool is_parallel;
};

}
#endif


