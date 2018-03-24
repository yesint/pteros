#include "pteros/analysis/task_base.h"
#include <thread>
#include "message_channel.h"
#include "pteros/core/pteros_error.h"
#include "data_container.h"
#include <iostream>

namespace pteros {

using Data_channel = Message_channel<std::shared_ptr<pteros::Data_container> > ;
using Data_channel_ptr = std::shared_ptr<Data_channel> ;

class Task_driver {
public:
    Task_driver(Task_base* _task);
    virtual ~Task_driver();
    void set_data_channel_and_system(const Data_channel_ptr& ch, const System &sys);
    void process_until_end();
    void process_until_end_in_thread ();
    void join_thread();        
    Frame_info get_last_info(){return last_info;}
private:
    Data_channel_ptr channel;
    Task_base* task;
    std::shared_ptr<Data_container> data;
    std::thread t;
    Frame_info last_info;
    bool stop_now; // Emergency stop flag for thread
    bool pre_process_done;
};


}
