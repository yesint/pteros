#ifndef TRAJ_FILE_READER_H
#define TRAJ_FILE_READER_H

#include "pteros/core/logging.h"
#include "pteros/analysis/options.h"
#include "message_channel.h"
#include "data_container.h"
#include <thread>

namespace pteros {

typedef Message_channel<std::shared_ptr<pteros::Data_container> > Data_channel;
typedef std::shared_ptr<Data_channel> Data_channel_ptr;


class Traj_file_reader {
public:
    Traj_file_reader(Options& options);

    bool is_frame_valid(int fr, float t);

    bool is_end_of_interval(int fr, float t);


    void run(const std::vector<std::string>& traj_files, const Data_channel_ptr& ch);

    ~Traj_file_reader();

    void join();

    void reader_thread_body(const std::vector<std::string>& traj_files, const Data_channel_ptr &channel);

private:
    int log_interval;
    float custom_start_time;
    float custom_dt;
    int first_frame, last_frame;
    float first_time, last_time;
    int skip;

    std::thread t;
    bool stop_now; // Emergency stop flag
    std::shared_ptr<spdlog::logger> log;
};

}

#endif // TRAJ_FILE_READER_H
