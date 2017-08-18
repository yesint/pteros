#include "pteros/analysis/task_base.h"
#include <thread>
#include "message_channel.h"
#include "pteros/core/pteros_error.h"
#include "data_container.h"
#include <iostream>

namespace pteros {

typedef Message_channel<std::shared_ptr<pteros::Data_container> > Data_channel;
typedef std::shared_ptr<Data_channel> Data_channel_ptr;

class Task_driver {
public:
    Task_driver(Task_base* _task);
    void set_data_channel(const Data_channel_ptr& ch);
    void init_with_first_frame(const System& sys);
    void process_first_frame();
    void process_until_end();
    void process_until_end_in_thread ();
    void join_thread();    
    void process_all(const System &sys);
    void process_all_in_thread(const System &sys);
    Frame_info get_last_info(){return last_info;}
private:
    Data_channel_ptr channel;
    Task_base* task;
    std::shared_ptr<Data_container> data;
    std::thread t;
    Frame_info last_info;
};


}
