#ifndef COMPILED_PLUGIN_BASE_H
#define COMPILED_PLUGIN_BASE_H

#include "pteros/analysis/consumer.h"
#include "pteros/analysis/trajectory_processor.h"
#include "boost/current_function.hpp"

namespace pteros {

struct Compiled_plugin_base: public Consumer {

    Compiled_plugin_base(Trajectory_processor* proc, const Options& opt): Consumer(proc){
        options = opt;
    }

    /// Unique textual label for this plugin instance
    /// Set by driver program
    std::string label;

    virtual std::string help(){
        return "\n\tThis plugin does not define any help information."
               "\n\tDig into the sources. Good luck :)";
    }

    void remove_jumps(Selection& sel, const Frame_info &info);

protected:
    /// Options for this particular plugin instance
    Options options;

    /// Handler overrides
    virtual void pre_process_handler(){
        try {
            Consumer::pre_process_handler();
        } catch(const Pteros_error& e){
            std::cout << std::endl << "(ERROR) Compiled plugin method: '"
                      << BOOST_CURRENT_FUNCTION << "'" << std::endl;
            std::cout << "(ERROR) in plugin instance: " << label << std::endl;
            std::cout << e.what() << std::endl;
            exit(1);
        }
    }

    virtual void post_process_handler(const Frame_info &info){
        try {
            Consumer::post_process_handler(info);
        } catch(const Pteros_error& e){
            std::cout << std::endl << "(ERROR) Compiled plugin method: '"
                      << BOOST_CURRENT_FUNCTION << "'" << std::endl;
            std::cout << "(ERROR) in plugin instance: " << label << std::endl;
            std::cout << e.what() << std::endl;
            exit(1);
        }
    }

    virtual void process_frame_handler(const Frame_info &info){
        try {
            Consumer::process_frame_handler(info);
        } catch(const Pteros_error& e){
            std::cout << std::endl << "(ERROR) Compiled plugin method: '"
                      << BOOST_CURRENT_FUNCTION << "'" << std::endl;
            std::cout << "(ERROR) Occured on frame " << info.valid_frame << std::endl;
            std::cout << "(ERROR) in plugin instance: " << label << std::endl;
            std::cout << e.what() << std::endl;
            exit(1);
        }
    }    

};

}

#endif
