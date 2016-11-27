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
protected:
    virtual bool is_parallel() = 0;
    virtual void pre_process() = 0;
    virtual void process_frame(const Frame_info& info) = 0;
    virtual void post_process(const Frame_info& info) = 0;


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
            return new _name(); \
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
