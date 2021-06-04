#include "pteros/analysis/task_base.h"
#include "task_driver.h"

using namespace std;
using namespace pteros;

TaskBase::TaskBase(): task_id(-1), n_consumed(0)
{
    //cout << "ctor: Task_base" << endl;
    driver.reset(new TaskDriver(this));
}

TaskBase::TaskBase(const TaskBase &other)
{
    //cout << "ctor copy: Task_base from " << other.task_id << endl;

    driver.reset(new TaskDriver(this));
    system = other.system;
    task_id = -1;
    n_consumed = 0;
}

void pteros::TaskBase::put_frame(const pteros::Frame &frame){
    system.frame(0) = frame;
}

void pteros::TaskBase::put_system(const pteros::System &sys){
    if(!system.num_atoms()) system = sys;
}
