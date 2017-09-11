#ifndef TASK_PLUGIN_H
#define TASK_PLUGIN_H

#include "pteros/analysis/task_base.h"
#include "pteros/analysis/options.h"
#include "pteros/analysis/jump_remover.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

namespace pteros {

/// Base class for plugins. Has options and jump remover
/// Also has constructor with options and set_options()

class Task_plugin: public Task_base {
public:

    Task_plugin(const Options& opt): options(opt), Task_base() {  }
    Task_plugin(const Task_plugin& other);

    Options options;
    Jump_remover jump_remover;

    Task_plugin* clone() const override {
        return nullptr;
    }

protected:

    void pre_process_handler() override
    {
        try {
            pre_process();
            jump_remover.remove_jumps(system); // Init jump remover

        } catch (const std::exception& e) {
            log->error("pre_process failed: {}", e.what());
            std::terminate();
        }
    }

    void process_frame_handler(const Frame_info& info) override
    {
        try {
            jump_remover.remove_jumps(system);
            process_frame(info);

        } catch (const std::exception& e) {
            log->error("process_frame failed on frame {}: {} ", info.valid_frame, e.what());
            std::terminate();
        }
    }

    virtual void post_process_handler(const Frame_info& info) override
    {
        try {
            post_process(info);

        } catch (const std::exception& e) {
            log->error("post_process failed: {} ", e.what());
            std::terminate();
        }
    }
};

}

#define _PLUGIN_(_name) \
    class _name: public Task_plugin { \
    public: \
        using Task_plugin::Task_plugin; \
        void set_id(int _id) override {\
            log = std::make_shared<spdlog::logger>(fmt::format(#_name ".{}",_id), Log::instance().console_sink); \
            log->set_pattern(Log::instance().generic_pattern); \
            task_id = _id; \
        } \

#define PLUGIN_PARALLEL(_name) \
    _PLUGIN_(_name) \
    _name* clone() const override { \
        return new _name(*this); \
    } \
    protected: \
        bool is_parallel() final { return true; }


#define PLUGIN_SERIAL(_name) \
    _PLUGIN_(_name) \
    protected: \
        bool is_parallel() final { return false; }



#endif // TASK_BASE_H
