#ifndef TASK_PLUGIN_H
#define TASK_PLUGIN_H

#include "pteros/analysis/task_base.h"
#include "pteros/analysis/options.h"
#include "pteros/analysis/jump_remover.h"

namespace pteros {

/// Base class for plugins. Has options and jump remover
/// Also has constructor with options and set_options()

class Task_plugin: public Task_base {
public:

    Task_plugin(const Options& opt): options(opt) { }    

    Options options;
    Jump_remover jump_remover;

protected:


    virtual void pre_process_handler(){
        pre_process();
        jump_remover.remove_jumps(system); // Init jump remover
    }


    virtual void process_frame_handler(const Frame_info& info){
        jump_remover.remove_jumps(system);
        process_frame(info);
    }
};

}


#define PLUGIN_PARALLEL(_name) \
    class _name: public Task_plugin { \
    public: \
        using Task_plugin::Task_plugin; \
        virtual _name* clone() const { \
            return new _name(*this); \
        } \
    protected: \
    virtual bool is_parallel() final { return true; }


#define PLUGIN_SERIAL(_name) \
    class _name: public Task_plugin { \
    public: \
        using Task_plugin::Task_plugin; \
        virtual _name* clone() const { \
            return new _name(*this); \
        } \
    protected: \
    virtual bool is_parallel() final { return false; }



#endif // TASK_BASE_H
