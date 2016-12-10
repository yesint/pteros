#ifndef TASK_BASE_H
#define TASK_BASE_H

#include "pteros/core/system.h"
#include "pteros/analysis/frame_info.h"

namespace pteros {

class Task_base {
    friend class Trajectory_reader;
public:
    Task_base(){}
    virtual ~Task_base(){}
    virtual Task_base* clone() const = 0;

    System system;
    int task_id;
protected:
    virtual bool is_parallel() = 0;
    virtual void pre_process() = 0;
    virtual void process_frame(const Frame_info& info) = 0;
    virtual void post_process(const Frame_info& info) = 0;

    // Handlers, which call actual functions
    // Could be overriden in subclasses
    virtual void pre_process_handler(){
        pre_process();
    }

    virtual void process_frame_handler(const Frame_info& info){
        process_frame(info);
    }

    virtual void post_process_handler(const Frame_info& info){
        post_process(info);
    }

private:

    void put_frame(const Frame& frame);
    void put_system(const System& sys);
};

}

#define TASK_DERIVED(_name) \
    class _name: public Task_base { \
    public: \
        _name(): Task_base() {} \
        virtual _name* clone() const { \
            return new _name(*this); \
        }

#define TASK_PARALLEL(_name) \
    TASK_DERIVED(_name) \
    protected: \
    virtual bool is_parallel() final { return true; }


#define TASK_SERIAL(_name) \
    TASK_DERIVED(_name) \
    protected: \
    virtual bool is_parallel() final { return false; }


#endif // TASK_BASE_H
