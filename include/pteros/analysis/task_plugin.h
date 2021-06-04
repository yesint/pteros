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

#include "pteros/analysis/task_base.h"
#include "pteros/analysis/options.h"
#include "pteros/analysis/jump_remover.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

namespace pteros {

/// Base class for plugins. Has options and jump remover
/// Also has constructor with options

class TaskPlugin: public TaskBase {
public:

    TaskPlugin(const Options& opt): options(opt), TaskBase() {  }
    TaskPlugin(const TaskPlugin& other);

    Options options;
    JumpRemover jump_remover;

    virtual std::string help(){ return ""; }

    TaskPlugin* clone() const override {
        return nullptr;
    }

protected:

    void before_spawn_handler() override;
    void pre_process_handler() override;
    void process_frame_handler(const FrameInfo& info) override;
    virtual void post_process_handler(const FrameInfo& info) override;
};

}

#define _TASK_(_name) \
    class _name: public TaskPlugin { \
    public: \
        using TaskPlugin::TaskPlugin; \
        void set_id(int _id) override {\
            log = create_logger(fmt::format(#_name ".{}",_id)); \
            task_id = _id; \
        } \

#define TASK_PARALLEL(_name) \
    _TASK_(_name) \
    _name* clone() const override { \
        return new _name(*this); \
    } \
    protected: \
        bool is_parallel() final { return true; }


#define TASK_SERIAL(_name) \
    _TASK_(_name) \
    protected: \
        bool is_parallel() final { return false; }

