#include "task_driver.h"

using namespace std;
using namespace pteros;

Task_driver::Task_driver(Task_base *_task): task(_task), stop_now(false)
{
    //cout << "ctor: Task_driver" << endl;
}

Task_driver::~Task_driver()
{
    if(t.joinable()){
        // Stop the thread
        stop_now = true;
        task->log->error("Ups! Stopping task driver thread on outer exception...");
        t.join();
    }
}

void Task_driver::set_data_channel_and_system(const Data_channel_ptr &ch, const System &sys){
    channel = ch;
    task->put_system(sys);
}

void Task_driver::process_until_end() {
    pre_process_done = false;
    while(channel->recieve(data)){
        if(stop_now) return; // Emergency stop point

        task->put_frame(data->frame);
        if(!pre_process_done){
            task->pre_process_handler();
            pre_process_done = true;
        }
        task->process_frame_handler(data->frame_info);
        ++task->n_consumed;
    }
    if(task->n_consumed>0){        
        task->post_process_handler(data->frame_info);        
    } else {
        task->log->warn("No frames consumed!");
    }
}

void Task_driver::process_until_end_in_thread() {
    t = std::thread(&Task_driver::process_until_end, this);
}

void Task_driver::join_thread() { t.join(); }

