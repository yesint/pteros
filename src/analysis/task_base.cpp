#include "pteros/analysis/task_base.h"
#include "task_driver.h"

using namespace std;
using namespace pteros;

Task_base::Task_base(): task_id(-1), n_consumed(0)
{
    //cout << "ctor: Task_base" << endl;
    driver.reset(new Task_driver(this));
}

Task_base::Task_base(const Task_base &other)
{
    //cout << "ctor copy: Task_base from " << other.task_id << endl;

    driver.reset(new Task_driver(this));
    system = other.system;
    task_id = -1;
    n_consumed = 0;
}

void pteros::Task_base::put_frame(const pteros::Frame &frame){
    system.frame(0) = frame;
}

void pteros::Task_base::put_system(const pteros::System &sys){
    if(!system.num_atoms()) system = sys;
}
