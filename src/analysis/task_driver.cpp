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

void Task_driver::set_data_channel(const Data_channel_ptr &ch){
    channel = ch;
}

void Task_driver::process_all(const System &sys) {
    init_with_first_frame(sys);
    process_first_frame();
    process_until_end();
}

void Task_driver::process_all_in_thread(const System &sys) {
    t = std::thread(&Task_driver::process_all, this, ref(sys));
}

void Task_driver::init_with_first_frame(const System &sys) {
    bool ok = channel->recieve(data);
    if(!ok) throw Pteros_error("Can't init instance of the task: no frames!");    
    task->put_system(sys);
    task->put_frame(data->frame);
    task->pre_process_handler();
}

void Task_driver::process_first_frame() {
    task->process_frame_handler(data->frame_info);
    ++task->n_consumed;
}

void Task_driver::process_until_end() {
    while(channel->recieve(data)){
        if(stop_now) return; // Emergency stop point

        task->put_frame(data->frame);
        task->process_frame_handler(data->frame_info);
        ++task->n_consumed;
    }
    if(task->n_consumed>0){
        // Only call post_process here for serial tasks!
        if(!task->is_parallel()) task->post_process_handler(data->frame_info);
        last_info = data->frame_info; // Save last processed frame
    } else {
        task->log->warn("No frames consumed!");
    }
}

void Task_driver::process_until_end_in_thread() {
    t = std::thread(&Task_driver::process_until_end, this);
}

void Task_driver::join_thread() { t.join(); }

