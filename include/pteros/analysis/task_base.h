#ifndef TASK_BASE_H
#define TASK_BASE_H

#include "pteros/core/system.h"
#include "pteros/analysis/frame_info.h"

// Forward declaration of the message channel
template<class T>
class Message_channel;

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

    virtual void set_id(int _id){ task_id = _id; }
    int get_id(){ return task_id; }

    System system;

    std::shared_ptr<spdlog::logger> log;

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

    int task_id;
    int n_consumed;

private:        

    void put_frame(const Frame& frame);
    void put_system(const System& sys);

    std::shared_ptr<Task_driver> driver;
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
