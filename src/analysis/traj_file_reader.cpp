#include "traj_file_reader.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/file_handler.h"
#include "pteros/core/utilities.h"

using namespace std;
using namespace pteros;

void process_suffix_value(const string& s, int* intval, float* floatval){
    size_t pos = s.find_last_of("0123456789");
    if(pos==string::npos) throw PterosError("A number with optional suffix required!");
    string val = s.substr(0,pos+1);
    string suffix = s.substr(pos+1);
    str_to_lower_in_place(val);
    str_to_lower_in_place(suffix);
    // Now analyze suffix
    if(intval!=nullptr && (suffix=="fr" || suffix=="")){
        *intval = str_to_int(val);
        if(floatval) *floatval = -1.0;
    } else if(floatval!=nullptr) {
        if(intval) *intval = -1;
        *floatval = str_to_float(val);
        if(suffix=="ps" || suffix=="t"){
            *floatval *= 1.0;
        } else if(suffix=="ns"){
            *floatval *= 1000.0;
        } else if(suffix=="us"){
            *floatval *= 1000000.0;
        } else if(suffix=="ms"){
            *floatval *= 1000000000.0;
        }
    } else if(floatval==nullptr && intval==nullptr){
        throw PterosError("Both int and float vals are NULL! WTF?");
    }
}

Traj_file_reader::Traj_file_reader(Options &options, int natoms){
    Natoms = natoms;

    // Separate reader logger (not registered since only used here)
    log = create_logger("traj_file");

    // Get parameters for reading frames
    process_suffix_value(options("b","-1").as_string(),
                         &first_frame, &first_time);
    process_suffix_value(options("e","-1").as_string(),
                         &last_frame, &last_time);
    // Skip interval
    skip = options("skip","-1").as_int();

    // Check for custom start time and dt
    process_suffix_value(options("t0","-1").as_string(),
                         nullptr, &custom_start_time);
    process_suffix_value(options("dt","-1").as_string(),
                         nullptr, &custom_dt);
    if(custom_start_time>=0 && custom_dt==-1) custom_dt = 1;
    if(custom_start_time==-1 && custom_dt>0) custom_start_time = 0;

    // Check if the range is valid
    if(first_frame>=0 && last_frame>=0 && last_frame<first_frame)
        throw PterosError("Last frame {} is smaller that first frame {}", last_frame,first_frame);
    if(first_time>=0 && last_time>=0 && last_time<first_time)
        throw PterosError("Last time {} is smaller that first time {}", last_time, first_time);

    log_interval = options("log","-1").as_int();
}

bool Traj_file_reader::is_frame_valid(int fr, float t){
    if(//------------------------------------------
            (first_frame<0 || fr>=first_frame)
            &&
            (first_time<0 || t>=first_time)
      ){//------------------------------------------
        // This frame is valid
        return true;
    } else {
        // This frame is invalid
        return false;
    }
}

bool Traj_file_reader::is_end_of_interval(int fr, float t){
    if(//------------------------------------------
            (fr>last_frame && last_frame>=0)
            ||
            (t>last_time && last_time>=0)
      ){//------------------------------------------
        // End reached
        return true;
    } else {
        return false;
    }
}

void Traj_file_reader::run(const vector<string> &traj_files, const DataChannel_ptr &ch){
    stop_now = false;
    t = std::thread( &Traj_file_reader::reader_thread_body, this, ref(traj_files), ref(ch) );
}

Traj_file_reader::~Traj_file_reader(){
    if(t.joinable()){
        // Stop the thread
        stop_now = true;
        log->error("Ups! Stopping reader thread on outer exception...");
        t.join();
    }
}

void Traj_file_reader::join(){ t.join(); }

void Traj_file_reader::reader_thread_body(const vector<string> &traj_files, const DataChannel_ptr &channel){
    try {
        int abs_frame = 0;
        float abs_time = 0.0;

        int valid_frame = -1;
        int frame_in_range = -1;
        // Saved first frame and time
        int first_valid_frame = -1;
        float first_valid_time = -1.0;

        bool finished = false;

        // Seek status:
        // 0 - don't need to seek
        // 1 - waiting for seeking
        int seek_status = 0;
        // Check if we need to seek for beginning
        if(first_frame>0 || first_time>0) seek_status = 1;

        for(const string& fname: traj_files){
            log->info("Reading trajectory {}...", fname);

            auto trj = FileHandler::open(fname,'r');

            // If we need to seek do it now if trajectory supports it
            if(seek_status==1 && trj->get_content_type().rand()){
                // Cast to random-access handler
                auto rand_trj = dynamic_cast<FileHandlerRandomAccess*>(trj.get());

                int last_fr;
                float last_t;
                rand_trj->tell_last_frame_and_time(last_fr,last_t);
                if(first_frame>0){
                    // If beyond this trajectory try the next one
                    if(first_frame>=last_fr){
                        log->info("First frame is {}, while this trajectory ends at {}.",first_frame,last_fr);
                        abs_frame += last_fr;
                        abs_time += last_t;
                        continue;
                    }
                    log->info("Fast forward to frame {}...",first_frame);
                    rand_trj->seek_frame(first_frame);
                } else if(first_time>0){
                    // If beyond this trajectory try the next one
                    if(first_time>=last_t){
                        log->info("First time is {}, while this trajectory ends at {}.",first_frame,last_fr);
                        abs_frame += last_fr;
                        abs_time += last_t;
                        continue;
                    }
                    log->info("Fast forward to time {}...",first_time);
                    rand_trj->seek_time(first_time);
                }
                seek_status = 0; // Seeking done
                // Set absolute frame count and time
                int fr;
                float t;
                rand_trj->tell_current_frame_and_time(fr,t);
                abs_frame += fr;
                abs_time += t;
                if(custom_dt>0) abs_time = custom_start_time + custom_dt*abs_frame;
            }

            --abs_frame;

            // Main loop over trajectory frames
            while(true){
                if(stop_now) return;

                // To avoid excessive copy operations we allocate a shared pointer
                // and will load data into its storage
                std::shared_ptr<DataContainer> data(new DataContainer);

                // Load data to this container
                bool good = trj->read(nullptr, &data->frame, FileContent().traj(true));

                // Check number of atoms
                if(data->frame.coord.size() != Natoms)
                    throw PterosError("Expected {} atoms but trajectory has {}.",data->frame.coord.size(),Natoms);

                // Check if EOF reached in trajectory
                if(!good) break;

                ++abs_frame; // Next absolute frame loaded

                // If time stamps are overriden, override time
                if(custom_dt>=0){
                    abs_time = custom_start_time + custom_dt*abs_frame;
                } else {
                    abs_time = data->frame.time;
                }

                if(log_interval>0 && abs_frame%log_interval==0)
                    log->info("At frame {}, {} ps",abs_frame,abs_time);

                // Check if end of requested interval is reached
                if( is_end_of_interval(abs_frame,abs_time) ){
                    // Send stop to the queue
                    channel->send_stop();
                    finished = true;
                    // exit loop
                    break;
                }

                // Check if new frame falls into needed range of time.
                // If not go to next frame
                if( !is_frame_valid(abs_frame,abs_time) ) continue;

                ++frame_in_range;

                // See if we need to skip it
                if(skip>0 && frame_in_range%skip!=0) continue;

                // This is valid frame
                ++valid_frame;

                if(valid_frame==0){
                    // This is the very first valid frame, set start time
                    first_valid_frame = abs_frame;
                    first_valid_time = abs_time;
                    // print a message
                    log->info("First valid frame is {}, {} ps",abs_frame,abs_time);
                }

                // Fill data container, which will be sent to the queue
                data->frame_info.absolute_time = abs_time;
                data->frame_info.absolute_frame = abs_frame;
                data->frame_info.valid_frame = valid_frame;
                data->frame_info.first_frame = first_valid_frame;
                data->frame_info.first_time = first_valid_time;
                data->frame_info.last_frame = abs_frame;
                data->frame_info.last_time = abs_time;

                // Send frame to the queue
                channel->send(data);

                // Do fast-forward skipping if asked
                if(skip>0){
                    if(trj->get_content_type().rand() && skip>0){
                        log->debug("Skipping {} frames by fast-forward...",skip);
                        try {
                            dynamic_cast<FileHandlerRandomAccess*>(trj.get())->seek_frame(abs_frame+skip);
                            abs_frame += skip;
                            frame_in_range += skip;
                        } catch(PterosError e){
                            log->debug("Can't seek, maybe EOF is reached");
                        }
                    }
                }
            } // Over frames

            log->info("Done with trajectory {}", fname);

            // If end reached break here too
            if(finished) break;

        } // Over trajectories

        // Send stop at the end
        channel->send_stop();

    } catch(const PterosError& e) {
        // Send stop if exception raised
        channel->send_stop();
        log->error(e.what());
    } catch(...) {
        log->critical("Some unknown terrible crash :-(");
        // Send stop if exception raised
        channel->send_stop();
    }
}
