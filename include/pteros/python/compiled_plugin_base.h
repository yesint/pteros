#ifndef COMPILED_PLUGIN_BASE_H
#define COMPILED_PLUGIN_BASE_H

#include "pteros/analysis/consumer.h"
#include "pteros/analysis/trajectory_processor.h"

namespace pteros {

struct Compiled_plugin_base: public Consumer {

    Compiled_plugin_base(Trajectory_processor* proc, Options_tree* opt): Consumer(proc){
        options = opt;
    }

    /// Unique textual label for this plugin instance
    /// Set by driver program
    std::string label;

    virtual string help(){
        return "\n\tThis plugin does not define any help information.\n\tDig into the sources. Good luck :)";
    }

protected:
    /// Options for this particular plugin instance
    Options_tree* options;

    /// Handler overrides
    virtual void pre_process_handler(){
        try {
            pre_process();
        } catch(Pteros_error e){
            cout << endl << "(ERROR) Compiled plugin method: '" << BOOST_CURRENT_FUNCTION << "'" << endl;
            cout << "(ERROR) in plugin instance: " << label << endl;
            e.print();
            exit(1);
        }
    }

    virtual void post_process_handler(const Frame_info &info){
        try {
            post_process(info);
        } catch(Pteros_error e){
            cout << endl << "(ERROR) Compiled plugin method: '" << BOOST_CURRENT_FUNCTION << "'" << endl;
            cout << "(ERROR) in plugin instance: " << label << endl;
            e.print();
            exit(1);
        }
    }

    virtual void process_frame_handler(const Frame_info &info){
        try {
            process_frame(info);
        } catch(Pteros_error e){
            cout << endl << "(ERROR) Compiled plugin method: '" << BOOST_CURRENT_FUNCTION << "'" << endl;
            cout << "(ERROR) Occured on frame " << info.valid_frame << endl;
            cout << "(ERROR) in plugin instance: " << label << endl;
            e.print();
            exit(1);
        }
    }

    virtual void window_started_handler(const Frame_info &info){
        try {
            window_started(info);
        } catch(Pteros_error e){
            cout << endl << "(ERROR) Compiled plugin method: '" << BOOST_CURRENT_FUNCTION << "'" << endl;
            cout << "(ERROR) Occured on frame " << info.valid_frame << endl;
            cout << "(ERROR) in plugin instance: " << label << endl;
            e.print();
            exit(1);
        }
    }


    virtual void window_finished_handler(const Frame_info &info){
        try {
            window_finished(info);
        } catch(Pteros_error e){
            cout << endl << "(ERROR) Compiled plugin method: '" << BOOST_CURRENT_FUNCTION << "'" << endl;
            cout << "(ERROR) Occured on frame " << info.valid_frame << endl;
            cout << "(ERROR) in plugin instance: " << label << endl;
            e.print();
            exit(1);
        }
    }

};

}

#endif
