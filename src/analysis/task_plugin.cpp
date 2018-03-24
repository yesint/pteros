#include "pteros/analysis/task_plugin.h"


pteros::Task_plugin::Task_plugin(const pteros::Task_plugin &other): Task_base(other)
{
    jump_remover = other.jump_remover;
}

void pteros::Task_plugin::before_spawn_handler() {
    before_spawn();
    // For parallel tasks init jump remover here
    if(is_parallel()) jump_remover.remove_jumps(system);
}

void pteros::Task_plugin::pre_process_handler()
{
    try {
        pre_process();

        // For serial tasks init jump remover here
        if(!is_parallel()) jump_remover.remove_jumps(system);

    } catch (const std::exception& e) {
        log->error("pre_process failed: {}", e.what());
        std::terminate();
    }
}

void pteros::Task_plugin::process_frame_handler(const pteros::Frame_info &info)
{
    try {
        jump_remover.remove_jumps(system);
        process_frame(info);

    } catch (const std::exception& e) {
        log->error("process_frame failed on frame {}: {} ", info.valid_frame, e.what());
        std::terminate();
    }
}

void pteros::Task_plugin::post_process_handler(const pteros::Frame_info &info)
{
    try {
        post_process(info);

    } catch (const std::exception& e) {
        log->error("post_process failed: {} ", e.what());
        std::terminate();
    }
}
