/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/


#ifndef OPTIONS_PARSER_H
#define OPTIONS_PARSER_H

#include <string>
#include <map>
#include <list>
#include <stack>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <pteros/core/pteros_error.h>
#include <json_spirit/json_spirit_writer_template.h>
#include <json_spirit/json_spirit_reader_template.h>

class Options_tree;

/// Single node
typedef boost::variant<
        boost::recursive_wrapper<Options_tree>,
        std::string,
        bool,
        int,
        double
    > Option_value;


class Options_tree {
    public:
        Options_tree(){
            value_iter = values.begin();
        }
        /// Construct empty object with given name
        Options_tree(std::string s){
            name = s;
            value_iter = values.begin();
        }
        /// Construct from json value
        Options_tree(json_spirit::mValue& json);

        /// Get all values of specified type from given key
        template<class T>
        std::list<T> get_values(std::string key);

        int count_options(std::string key);

        /// Get single value. Only one value should exist for success.
        template<class T>
        T get_value(std::string key);

        /// Get single value. Only one value should exist for success.
        /// If no value found return default value
        template<class T>
        T get_value(std::string key, T def_val);

        template<class T>
        void add_value(std::string key, T val);

        /// Get single option (only one option with given key should exist)
        Options_tree& get_option(std::string key);

        /// Get all options with given key
        std::list<Options_tree*> get_options(std::string key);

        /// Convert to json
        json_spirit::mValue to_json();
        std::string to_json_string();

        /// Read from json
        void from_json(json_spirit::mValue& json);

        /** Create tree from the command line
            the syntax of nested options is:
            --name1 val1 val2 val3 --name2 v1 v2 --end-name2 val 4 --end-name1
            Options without trailing --end become children of the previous
            option with --end:
            --name1 v1 --nested1 1 --nested2 2 --end-name1
            If no --ends at all then all options are children of root node.
            Root node is added automatically if needed.
        */
        void from_command_line(std::string cmd);
        void from_command_line(std::vector<std::string>& tokens);
        void from_command_line(int argc, char** argv);
        /// Convert to command line representation
        std::string to_command_line();

        /// Create options tree from idented input
        void from_indented(std::string str);
        std::string to_indented(int level = 0);

        /// >> extracts next value from given option
        template <class T>
        Options_tree& operator>>(T& val);
        /// Specializations for int and double
        Options_tree& operator>>(int& val);
        Options_tree& operator>>(double& val);

    protected:

        std::string name; //Name of option
        std::list<Option_value> values; //Values of option
        // Finds given key
        void find_key(std::string& key, std::list<Options_tree*>& matches, bool create_new = false);

        // Find single options
        Options_tree* get_single_option(std::string& key);

        // Recursively transforms json to Options_tree
        void json_to_tree(json_spirit::mValue& json, Options_tree& tree);
        // Itartor of values for >> operator
        std::list<Option_value>::iterator value_iter;
};

using namespace std;

template<class T>
T Options_tree::get_value(std::string key){    
    // Find matching options
    list<Options_tree*> v;
    find_key(key,v);
    // If more then one Options_tree matches, throw error
    if(v.size()>1) throw pteros::Pteros_error("More then one Options_tree node matches key '"+key+"'!");
    // If no such node also throw an error
    if(v.size()==0) throw pteros::Pteros_error("No Options_tree nodes match key '"+key+"'!");
    // We need exactly one values of type T
    T* ptr = NULL;
    T* cur;
    int k = 0;
    for(Option_value& val : v.front()->values){
         cur = boost::get<T>(&val);
        if(cur){
            ++k;
            ptr = cur;
        }
        if(k>1) throw pteros::Pteros_error("Exactly one value of requested type should be present "
                                           " for key '" + key + "'!");
    }
    if(!ptr) throw pteros::Pteros_error("No values of requested type are found in key '"+key+"'!");
    return *ptr;
}

template<class T>
T Options_tree::get_value(std::string key, T def_val){    
    // Find matching options
    list<Options_tree*> v;
    find_key(key,v);
    // If more then one Options_tree matches, throw error
    if(v.size()>1) throw pteros::Pteros_error("More then one Options_tree matches key '"+key+"'!");
    if(v.size()==0) return def_val;
    // We need exactly one values of type T
    T* ptr = NULL;
    T* cur;
    int k = 0;
    for(Option_value& val: v.front()->values){
        cur = boost::get<T>(&val);
        if(cur){
            ++k;
            ptr = cur;
        }
        if(k>1) throw pteros::Pteros_error("Exactly one value of requested type should be present "
                                           " for key '" + key + "'!");
    }
    if(ptr){
        return *ptr;
    } else
        return def_val;
}

template<class T>
void Options_tree::add_value(std::string key, T val){

    // Find matching options
    list<Options_tree*> v;
    find_key(key,v, true);
    //cout << boost::get<Options_tree*>(v.front())->name << endl;
    // If more then one Options_tree matches, throw error
    if(v.size()>1) throw pteros::Pteros_error("More then one Options_tree matches key '"+key+"'!");

    // Now add value
    v.front()->values.push_back(val);
}


template<class T>
list<T> Options_tree::get_values(std::string key){
    list<T> res;
    // Find matching options
    list<Options_tree*> v;
    find_key(key,v);
    // If more then one Options_tree matches, throw error
    if(v.size()>1) throw pteros::Pteros_error("More then one Options_tree matches!");
    // Return all values, which are of type T
    T* ptr;
    for(Option_value& val: v.front()->values){
        ptr = boost::get<T>(&val);
        if(ptr) res.push_back(*ptr);
    }
    return res;
}

template<class T>
bool add_to_json_array(json_spirit::mArray& arr, Option_value& o){
    T* ptr = boost::get<T>(&o);
    if(ptr){
        arr.push_back(*ptr);
        return true;
    } else
        return false;
}

template<class T>
bool add_to_command_line(std::string& cmd, Option_value& o){
    T* ptr = boost::get<T>(&o);
    if(ptr){
        cmd += " " + boost::lexical_cast<string>(*ptr);
        return true;
    } else
        return false;
}

template<class T>
bool add_to_indented(std::string& str, Option_value& o, int level){
    T* ptr = boost::get<T>(&o);
    if(ptr){
        str += string(level,' ') + boost::lexical_cast<string>(*ptr) + '\n';
        return true;
    } else
        return false;
}


template <class T>
Options_tree& Options_tree::operator>>(T& val){
    if(value_iter==values.end())
        throw pteros::Pteros_error("No more values!");
    // Set val to current value
    try {
        val = boost::get<T>(*value_iter);
    } catch(boost::bad_get e){
        cout << e.what() << endl;
    }

    // Advance iterator
    value_iter++;
    return *this;
}



#endif
