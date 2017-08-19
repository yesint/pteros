#include "pteros/analysis/task_plugin.h"


pteros::Task_plugin::Task_plugin(const pteros::Task_plugin &other): Task_base(other)
{
    jump_remover = other.jump_remover;
}
