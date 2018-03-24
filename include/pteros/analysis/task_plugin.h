#ifndef TASK_PLUGIN_H
#define TASK_PLUGIN_H

#include "pteros/analysis/task_base.h"
#include "pteros/analysis/options.h"
#include "pteros/analysis/jump_remover.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

namespace pteros {

/// Base class for plugins. Has options and jump remover
/// Also has constructor with options

class Task_plugin: public Task_base {
public:

    Task_plugin(const Options& opt): options(opt), Task_base() {  }
    Task_plugin(const Task_plugin& other);

    Options options;
    Jump_remover jump_remover;

    virtual std::string help(){ return ""; }

    Task_plugin* clone() const override {
        return nullptr;
    }

protected:

    void before_spawn_handler() override;
    void pre_process_handler() override;
    void process_frame_handler(const Frame_info& info) override;
    virtual void post_process_handler(const Frame_info& info) override;
};

}

#define _TASK_(_name) \
    class _name: public Task_plugin { \
    public: \
        using Task_plugin::Task_plugin; \
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



#endif // TASK_BASE_H
