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
template<class T> class Message_channel;

namespace pteros {

// Forward declarations
class Data_container;
class Task_driver;


class Task_base {
    friend class Task_driver;
    friend class Trajectory_reader;

public:
    Task_base();
    Task_base(const Task_base& other);
    virtual ~Task_base(){}
    virtual Task_base* clone() const = 0;

    int get_id(){ return task_id; }

    System system;

    std::shared_ptr<spdlog::logger> log;

    virtual void pre_process() = 0;
    virtual void process_frame(const Frame_info& info) = 0;
    virtual void post_process(const Frame_info& info) = 0;

    // Default implementation of collector for parallel tasks
    virtual void collect_data(const std::vector<std::shared_ptr<Task_base>>& tasks, int n_frames){}

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

    virtual void process_frame_handler(const Frame_info& info){
        process_frame(info);
    }

    virtual void post_process_handler(const Frame_info& info){
        post_process(info);
    }

    int task_id;
    int n_consumed;

private:        

    void put_frame(const Frame& frame);
    void put_system(const System& sys);

    std::shared_ptr<Task_driver> driver;
};

}
