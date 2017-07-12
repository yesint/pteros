#ifndef TASK_BASE_H
#define TASK_BASE_H

#include "pteros/core/system.h"
#include "pteros/analysis/frame_info.h"

// Forward declaration of the message channel
template<class T>
class Message_channel;


namespace pteros {

// Forward declaration
class Data_container;

typedef Message_channel<std::shared_ptr<pteros::Data_container> > Data_channel;
typedef std::shared_ptr<Data_channel> Data_channel_ptr;


class Task_base {
    friend class Task_driver;
    friend class Trajectory_reader;
public:
    Task_base(){}
    virtual ~Task_base(){}
    virtual Task_base* clone() const = 0;

    System system;
    int task_id;
    int n_consumed;
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
    Data_channel_ptr channel;

    void put_frame(const Frame& frame);
    void put_system(const System& sys);

    void set_channel(Data_channel_ptr ch);

    void init_with_first_frame();
    void process_first_frame();
    void consume_until_end();
    void consume_until_end_in_thread();
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
