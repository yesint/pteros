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


#pragma once

#include "pteros/core/system.h"
#include "pteros/analysis/frame_info.h"
#include <spdlog/spdlog.h>

// Forward declaration of the message channel
template<class T> class MessageChannel;

namespace pteros {

// Forward declarations
class DataContainer;
class TaskDriver;


class TaskBase {
    friend class TaskDriver;
    friend class TrajectoryReader;

public:
    TaskBase();
    TaskBase(const TaskBase& other);
    virtual ~TaskBase(){}
    virtual TaskBase* clone() const = 0;

    int get_id(){ return task_id; }

    System system;

    std::shared_ptr<spdlog::logger> log;

    virtual void pre_process() = 0;
    virtual void process_frame(const FrameInfo& info) = 0;
    virtual void post_process(const FrameInfo& info) = 0;

    // Default implementation of collector for parallel tasks
    virtual void collect_data(const std::vector<std::shared_ptr<TaskBase>>& tasks, int n_frames){}

    // Default implementation of global preprocess for parallel tasks
    virtual void before_spawn(){}

protected:
    virtual void set_id(int _id){ task_id = _id; }

    virtual bool is_parallel() = 0;    

    // Handlers, which call actual functions
    // Could be overriden in subclasses
    virtual void before_spawn_handler(){
        before_spawn();
    }

    virtual void pre_process_handler(){
        pre_process();
    }

    virtual void process_frame_handler(const FrameInfo& info){
        process_frame(info);
    }

    virtual void post_process_handler(const FrameInfo& info){
        post_process(info);
    }

    int task_id;
    int n_consumed;

private:        

    void put_frame(const Frame& frame);
    void put_system(const System& sys);

    std::shared_ptr<TaskDriver> driver;
};

}
