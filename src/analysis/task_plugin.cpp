#include "pteros/analysis/task_plugin.h"


pteros::TaskPlugin::TaskPlugin(const pteros::TaskPlugin &other): TaskBase(other)
{
    jump_remover = other.jump_remover;
}

void pteros::TaskPlugin::before_spawn_handler() {
    before_spawn();
    // For parallel tasks init jump remover here
    if(is_parallel()) jump_remover.remove_jumps(system);
}

void pteros::TaskPlugin::pre_process_handler()
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

void pteros::TaskPlugin::process_frame_handler(const pteros::FrameInfo &info)
{
    try {
        jump_remover.remove_jumps(system);
        process_frame(info);

    } catch (const std::exception& e) {
        log->error("process_frame failed on frame {}: {} ", info.valid_frame, e.what());
        std::terminate();
    }
}

void pteros::TaskPlugin::post_process_handler(const pteros::FrameInfo &info)
{
    try {
        post_process(info);

    } catch (const std::exception& e) {
        log->error("post_process failed: {} ", e.what());
        std::terminate();
    }
}
