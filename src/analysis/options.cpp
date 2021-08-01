/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#include "pteros/analysis/options.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/utilities.h"
#include <charconv>
#include <sstream>

using namespace std;
using namespace pteros;

namespace pteros {

void parse_command_line(int argc, char** argv,
                        Options& toplevel,
                        std::string task_tag,
                        std::vector<Options>& tasks){

    if(argc<2) return; // Command line is empty

    bool in_task = false; // start with toplevel
    bool has_key = false;
    // tmp option
    Option o;
    o.name = "";
    // tmp task
    Options tsk;
    tsk.task_name = "";

    for(int i=1;i<argc;++i){        
        string str(argv[i]);

        if(str=="-") throw PterosError("Lone '-'!");

        //See if token contains '-' as first symbol not followed by number
        if(str[0]=='-' && str.size()>1 && !isdigit(str[1]) ){
            // This is key
            has_key = true;
            // strip '-'
            str = str.substr(1);

            // If we have filled option already, add it
            if(o.name!=""){
                if(!in_task){
                    toplevel.data.push_back(o);
                } else {
                    tsk.data.push_back(o);
                }
            }

            // See if we got task tag
            if(str==task_tag){
                if(i==argc-1) throw PterosError("Incomplete task at the end of command line!");
                // If we already been in task mode, then end old task
                if(in_task){
                    if(tsk.task_name=="") throw PterosError("Task without name after {}!",argv[i-1]);
                    tasks.push_back(tsk);
                }
                in_task = true;
                has_key = false;
                o.name = "";
                o.data.clear();
                tsk.task_name = "";
                tsk.data.clear();
                continue;
            }

            // set name
            o.name = str;
            // clear values
            o.data.clear();
        } else {
            // This is value
            if(!has_key && !in_task) throw PterosError("Error! Found value '{}' without option name!",str);
            if(!has_key && in_task){
                // This is task name
                if(tsk.task_name!="") throw PterosError("Error: double task name '{}'",str);
                tsk.task_name = str;
                if(i==argc-1) tasks.push_back(tsk);
                continue;
            }
            // Add it
            o.data.push_back(str);
        }
    }

    //TODO: last task without options doesn't work...
    // Hm, not sure if this is still the case...

    // At the end see where to put last option
    if(o.name!=""){
        if(!in_task){
            toplevel.data.push_back(o);
        } else {
            tsk.data.push_back(o);
            // Add task itself
            tasks.push_back(tsk);
        }
    }    
}

// without nested tasks
void parse_command_line(int argc, char** argv, Options& toplevel){
    vector<Options> dum;
    parse_command_line(argc,argv,toplevel,"",dum);
}

} // namespece

#ifdef DEBUG
void Options::debug(){
    cout << "task_name = " << task_name << endl;
    for(auto& o: data){
        o.debug();
    }
}

void Option::debug(){
    cout << name << ": ";
    for(auto& o: data){
        cout << o << " ";
    }
    cout << endl;
}
#endif

const Option& Options::operator()(std::string key) const {
    int ind = -1;
    int found = 0;
    for(int i=0;i<data.size();++i){
        if(data[i].name==key){
            ind = i;
            ++found;
        }
    }
    if(found==0) throw PterosError("Key '{}' not found!", key);
    if(found>1) throw PterosError("More than one key '{}' found!", key);
    return data[ind];
}

const Option& Options::operator()(std::string key, std::string default_val) {
    int ind = -1;
    int found = 0;
    for(int i=0;i<data.size();++i){
        if(data[i].name==key){
            ind = i;
            ++found;
        }
    }
    if(found>1) throw PterosError("More than one key '{}' found!",key);
    if(found==0 || data[ind].data.empty()){
        // Default value case
        Option tmp;
        tmp.name = key;        
        // Split default_val by " "
        stringstream ss(default_val);
        string buf;
        while(ss >> buf) tmp.data.push_back(buf);        
        // If tmp.data is still empty here, than default_val is empty string or a space(s)
        // In this case we add it as is
        if(tmp.data.empty()) tmp.data.push_back(default_val);

        // Default value is added to data and returned as usually
        data.push_back(tmp);
        return data.back();
    }
    // Normal case
    return data[ind];
}

bool pteros::Options::has(string key)
{
    int ind = -1;
    int found = 0;
    for(int i=0;i<data.size();++i){
        if(data[i].name==key){
            ind = i;
            ++found;
        }
    }
    if(found>1) throw PterosError("More than one key '{}' found!",key);
    if(found==0)
        return false;
    else
        return true;
}


//----------------------------------------------------------------

int Option::as_int() const {
    if(data.size()!=1) throw PterosError("Only one INT value expected for key '{}'!",name);
    return str_to_int(data[0]);
}

float Option::as_float() const {
    if(data.size()!=1) throw PterosError("Only one FLOAT value expected for key '{}'!",name);
    return str_to_float(data[0]);
}

bool Option::as_bool() const {
    if(data.size()!=1) throw PterosError("Only one BOOL value expected for key '{}'!",name);
    string s(data[0]);
    str_to_lower_in_place(s);

    if(s=="true"){
        return true;
    } else if(s=="false"){
        return false;
    } else {
        throw PterosError("Value '{}' is not BOOL!",s);
    }

}

string Option::as_string() const {        
    if(data.size()!=1) throw PterosError("Only one STRING value expected for key '{}'!",name);
    return data[0];
}

vector<int> Option::as_ints() const {
    if(data.empty()) throw PterosError("One or more INT values are expected for key '{}'!",name);
    vector<int> res;
    for(auto& str: data){
        res.push_back( str_to_int(str) );
    }
    return res;
}

vector<float> Option::as_floats() const {
    if(data.empty()) throw PterosError("One or more FLOAT values are expected for key '{}'!",name);
    vector<float> res;
    for(auto& str: data){
        res.push_back( str_to_float(str) );
    }
    return res;
}

Eigen::VectorXf Option::as_VectorXf() const {
    if(data.empty()) throw PterosError("One or more FLOAT values are expected for key '{}'!",name);
    vector<float> vec;
    for(auto& str: data){
        vec.push_back( str_to_float(str) );
    }
    Eigen::VectorXf res(vec.size());
    for(int i=0;i<vec.size();++i) res(i)=vec[i];
    return res;
}


vector<bool> Option::as_bools() const {
    if(data.empty()) throw PterosError("One or more BOOL values are expected for key '{}'!",name);
    vector<bool> res;
    for(auto s: data){
        str_to_lower_in_place(s);
        if(s=="true"){
            res.push_back(true);
        } else if(s=="false"){
            res.push_back(false);
        } else {
            throw PterosError("Value '{}' is not BOOL!",s);
        }
    }
    return res;
}

vector<string> Option::as_strings() const {
    if(data.empty()) throw PterosError("One or more STRING values are expected for key '{}'!",name);
    return data;
}




