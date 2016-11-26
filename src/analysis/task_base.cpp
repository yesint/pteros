#include "pteros/analysis/task_base.h"


void pteros::Task_base::put_frame(const pteros::Frame &frame){
    system.Frame_data(0) = frame;
}

void pteros::Task_base::put_system(const pteros::System &sys){
    system = sys;
}
