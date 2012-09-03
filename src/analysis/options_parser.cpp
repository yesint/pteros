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


#include "pteros/analysis/options_parser.h"
#include <fstream>

using namespace std;

bool is_object(Option_value& o){
    Options_tree* ptr = boost::get<Options_tree>(&o);
    if(ptr)
        return true;
    else
        return false;
}


Options_tree* Options_tree::get_single_option(std::string& key){
    // Find matching options
    list<Options_tree*> v;
    find_key(key,v);
    // If more then one Options_tree matches, throw error
    if(v.size()>1) throw pteros::Pteros_error("More then one options '"+key+" are found'!");
    if(v.size()==0) throw pteros::Pteros_error("Can't find required option '"+key+"'!");
    return v.front();
}

std::list<Options_tree*> Options_tree::get_options(std::string key){
    // Find matching options
    list<Options_tree*> v;
    find_key(key,v);
    return v;
}

int Options_tree::count_options(std::string key){
    // Find matching options
    list<Options_tree*> v;
    find_key(key,v);
    return v.size();
}


void Options_tree::find_key(std::string& key,
                            list<Options_tree*>& matches, bool create_new){
    Options_tree *ptr;
    // Split key by / and add tokens to path
    list<string> path;
    boost::split(path,key,boost::is_any_of("./"),boost::token_compress_on);

    matches.clear(); // No matche so far
    matches.push_back(this); // Include current key to search

    // If empty key is given return current node, which was just added
    if(path.size()==0) return;

    list<Options_tree*> new_matches;

    // Cycle over elements in path
    BOOST_FOREACH(string& s, path){
        // Find name s in the list of nodes given by matches

        // Ignore empty names
        if(s=="") continue;

        if(matches.size()){
            // For each element in matches search for nodes child with name s
            // All found nodes are returned as new_mathes
            new_matches.clear();
            BOOST_FOREACH(Options_tree* cur, matches){
                // Cycle over values of cur
                BOOST_FOREACH(Option_value& val, cur->values){
                    ptr = boost::get<Options_tree>(&val);
                    // If not NULL then it is an Options_tree and we can check name
                    if(ptr){
                        if(ptr->name == s){
                            // Name matches
                            new_matches.push_back(ptr);
                        }
                    }
                }
            }

            // If asked to crate node, create node names s
            if(create_new && new_matches.size()==0){
                if(matches.size()>1) throw "Can't create child node "+s+"!";
                Options_tree o;
                o.name = s;
                matches.front()->values.push_back(o);
                new_matches.push_back( boost::get<Options_tree>(&matches.front()->values.back()) );
            }

            // Reset matches to new found list nodes
            matches = new_matches;

        } else {            
            break;
        }

    }
    // Set value_iter for all matches
    BOOST_FOREACH(Options_tree* o,matches){
        o->value_iter = o->values.begin();
    }        
}

string Options_tree::to_json_string(){
    json_spirit::mValue val = to_json();
    return json_spirit::write_string(val,true);
}

json_spirit::mValue Options_tree::to_json(){
    using namespace json_spirit;

    mArray vals;
    BOOST_FOREACH(Option_value& o,values){
        {// Try to get nested Options_tree
            Options_tree* ptr = boost::get<Options_tree>(&o);
            if(ptr){
                mObject obj;
                obj[ptr->name] = ptr->to_json();
                vals.push_back(obj);
                continue;
            }
        }

        if( add_to_json_array<double>(vals,o) )
            continue;
        if( add_to_json_array<int>(vals,o) )
            continue;
        if( add_to_json_array<string>(vals,o) )
            continue;
        if( add_to_json_array<bool>(vals,o) )
            continue;
        throw pteros::Pteros_error("Can't recognize option type");
    }

    return mValue(vals);
}

void Options_tree::json_to_tree(json_spirit::mValue& json, Options_tree& tree){
    using namespace json_spirit;

    if(json.type() == obj_type){
        // If this is object cycle over fields
        pair<string,mValue> e;
        BOOST_FOREACH(e, json.get_obj()){
            // Create new tree node with name of this field
            Options_tree t(e.first);
            // Recurse for value
            json_to_tree(e.second,t);
            // Add this node as child
            tree.values.push_back(t);
        }
    } else if(json.type() == array_type){
        // If this is array simply recurse for each element
        BOOST_FOREACH(mValue& el, json.get_array()){
            json_to_tree(el,tree);
        }
    } else if(json.type() == int_type){
        tree.values.push_back(json.get_int());
    } else if(json.type() == real_type){
        tree.values.push_back(json.get_real());
    } else if(json.type() == str_type){
        tree.values.push_back(json.get_str());
    } else if(json.type() == bool_type){
        tree.values.push_back(json.get_bool());
    }
}

Options_tree::Options_tree(json_spirit::mValue& json){
    // Read json
    name = "root";
    values.clear();
    json_to_tree(json,*this);
}

void Options_tree::from_json(json_spirit::mValue& json){
    name = "root";
    values.clear();
    json_to_tree(json,*this);
}


void Options_tree::from_command_line(std::vector<std::string>& tokens){
    // If no tokens, just exit
    if(tokens.empty()) return;

    // If first token is special option --json, then treat the second token as
    // a json string. The second token should be in "" to be treated as single token
    if(tokens[0]=="--json"){
        json_spirit::mValue v;
        // Check if we can open tokens[1] as file
        ifstream f(tokens[1].c_str());
        if(f){
            cout << "Reading JSON from file '" << tokens[1] << "'..." << endl;
            json_spirit::read_stream(f,v);
        } else {
            cout << "Reading JSON from command line..."<<endl;
            stringstream s(tokens[1]);
            json_spirit::read_stream(s,v);
        }
        from_json(v);
        return;
    }


    // Create root node
    name = "root";
    values.clear();

    stack<Options_tree*> cur; //Stack of nested options

    // Add root
    cur.push(this);

    bool is_closed;

    // Parse input
    for(int i=0; i<tokens.size(); ++i){
        string tok = tokens[i]; // Current token

        // End of group
        if(tok == "]"){
            cur.pop();
            continue;
        }

        // Beginning of the group
        if(tok == "["){
            // See if it is a group by searching for [ in the next token
            int level = 1;
            // This is a group, so seacrh for matching closing ]
            for(int j=i+1; j<tokens.size(); ++j){
                if(tokens[j]=="[") level++;
                if(tokens[j]=="]"){
                    level--;
                    if(level<0) throw pteros::Pteros_error("Unbalanced ']' found!");
                    if(level==0) break; // Closed found
                }
            }
            if(level) throw pteros::Pteros_error("Unbalanced '[' found!");

            // Pointer to parent object
            Options_tree* ptr = boost::get<Options_tree*>(cur.top());
            cur.push(ptr); // Push it to cur again
            // All subsequent values are added to parent object
            continue;
        }

        // If token starts with -- then this is a name
        if(tok.substr(0,2)=="--"){
            // This is begin of new option
            // pop last option if this is not root
            if(cur.size()>1) cur.pop();

            Options_tree opt;
            opt.name = tok.substr(2);
            cur.top()->values.push_back(opt);

            //cout << "==> " << opt.name << " added to "
            //    << boost::get<Options_tree*>(cur.top())->name << endl;

            // Pointer to last added object
            Options_tree* ptr = boost::get<Options_tree>(&cur.top()->values.back());
            cur.push(ptr); // Push it to cur
            // All subsequent values are added to new object

        } else {
            // This is value
            // Convert types here!
            if(tok=="true"){
                cur.top()->values.push_back(true);
            } else if(tok=="false") {
                cur.top()->values.push_back(false);
            } else {
                try{
                    int i = boost::lexical_cast<int>(tok);
                    cur.top()->values.push_back(i);
                } catch(boost::bad_lexical_cast) {
                    try {
                        double i = boost::lexical_cast<double>(tok);
                        cur.top()->values.push_back(i);
                    } catch(boost::bad_lexical_cast) {
                        cur.top()->values.push_back(tok);
                    }
                }
            } // Type conversion

            //cout << "\t==> " << tok << " added to "
            //    << boost::get<Options_tree*>(cur.top())->name << endl;

        }
    }

    // If our top root node has only one child, which is an option
    // then substitute top root
    /*
    if(values.size()==1 && ::is_object(values.front())){
        Options_tree p = boost::get<Options_tree>(values.front());
        *this = p;
    }
    */

}

void Options_tree::from_command_line(std::string cmd){
    // Separate [ and ] by spaces
    boost::replace_all(cmd,"["," [ ");
    boost::replace_all(cmd,"]"," ] ");
    vector<string> tokens;
    // Split a string into the tokens
    boost::escaped_list_separator<char> sep('\\',' ','\"');
    boost::tokenizer<boost::escaped_list_separator<char> > tok(cmd,sep);
    boost::tokenizer<boost::escaped_list_separator<char> >::iterator t;
    for(t=tok.begin();t!=tok.end();++t){
        if(*t!="") tokens.push_back(*t);
    }

    from_command_line(tokens);
}

void Options_tree::from_command_line(int argc, char** argv){
    std::vector<std::string> tok;
    std::vector<std::string> parts;
    for(int i=1;i<argc;++i){ // We don't need name of program
        // Parse [ and ]
        string s = "";
        BOOST_FOREACH(char ch, argv[i]){
            if(ch!='[' && ch!=']'){
                s += ch;
            } else {
                if(s!="") tok.push_back(s);
                s = ch;
                tok.push_back(s);
                s = "";
            }
        }
        if(s!="") tok.push_back(s);

        //tok.push_back(string(argv[i]));
    }

    //for(int i=0;i<tok.size();++i) cout << tok[i] << endl;


    from_command_line(tok);
}

string Options_tree::to_command_line(){
    string cmd;

    cmd += "--"+name;

    bool need_end = false;
    bool has_nested = false;

    BOOST_FOREACH(Option_value& o,values){
        {// Try to get nested Options_tree
            Options_tree* ptr = boost::get<Options_tree>(&o);
            if(ptr){
                cmd+= " [ "+ptr->to_command_line();
                has_nested = true;
                continue;
            }
        }

        if( add_to_command_line<int>(cmd,o) )
            continue;
        if( add_to_command_line<double>(cmd,o) )
            continue;
        if( add_to_command_line<string>(cmd,o) )
            continue;
        if( add_to_command_line<bool>(cmd,o) )
            continue;
    }

    if(has_nested) cmd += " ] ";
    return cmd;
}

void Options_tree::from_indented(string str){
    // Indentation level
    int level = 0;
    int i;
    //First split the input by \n
    vector<string> lines;
    boost::split(lines,str,boost::is_any_of("\n\r"),boost::token_compress_on);

    // Create root node
    name = "root";
    values.clear();

    stack<Options_tree*> cur; //Stack of nested options
    stack<int> levels; // Stack of associated levels

    cur.push(this); // Add root
    levels.push(-1);

    Options_tree* ptr;
    int int_val;
    float float_val;

    bool added = false;

    BOOST_FOREACH(string& line, lines){
        // Analyse this line.
        // First determine the indentation
        i = 0;
        while(line[i]==' ') ++i;
        // If this is already end of line skip
        if(i==line.size()) continue;

        if(added){
            added = false;
            levels.push(i);
            level = i;
        } else if(i<level) {
            while(i<levels.top()){
                cur.pop();
                levels.pop();
            }
            level = i;
        }

        //cout << line << " indent: " << i << " level: " << level << " current: " << cur.top()->name << endl;

        // Add option to current root
        // Determine what we have
        if(line[i]=='"'){
            // This is string value
            cur.top()->values.push_back(line.substr(i+1,line.find_first_of('"',i+1)-i-1));
        } else if(line.substr(i)=="true") {
            // true
            cur.top()->values.push_back(true);
        } else if(line.substr(i)=="false") {
            // false
            cur.top()->values.push_back(false);
        } else {
            try {
                float_val = boost::lexical_cast<float>(line.substr(i));
                // This is float value
                cur.top()->values.push_back(float_val);
            } catch(boost::bad_lexical_cast) {
                try {
                    int_val = boost::lexical_cast<int>(line.substr(i));
                    // This is int value
                    cur.top()->values.push_back(int_val);
                } catch(boost::bad_lexical_cast) {
                    // If we are here this is an option name
                    Options_tree o(line.substr(i));
                    cur.top()->values.push_back(o);

                    ptr = boost::get<Options_tree>(&cur.top()->values.back());
                    cur.push(ptr); // Push it to cur

                    //cout << "Added option "+cur.top()->name << endl;

                    // Indent is now increased
                    added = true;
                }
            }
        }

    }
}

template<>
bool add_to_indented<string>(std::string& str, Option_value& o, int level){
    string* ptr = boost::get<string>(&o);
    if(ptr){
        str += string(level,' ') + '"' + (*ptr) + "\"\n";
        return true;
    } else
        return false;
}

string Options_tree::to_indented(int level){
    string str;

    if(name!="root"){
        str += string(level,' ') + name + "\n";
        ++level;
    }

    BOOST_FOREACH(Option_value& o,values){
        {// Try to get nested Options_tree
            Options_tree* ptr = boost::get<Options_tree>(&o);
            if(ptr){
                str += ptr->to_indented(level);
                continue;
            }
        }

        if( add_to_indented<int>(str,o,level) )
            continue;
        if( add_to_indented<double>(str,o,level) )
            continue;
        if( add_to_indented<string>(str,o,level) )
            continue;
        if( add_to_indented<bool>(str,o,level) )
            continue;
    }

    return str;
}


Options_tree& Options_tree::get_option(std::string key){
    // Set iterator to the first value to reset >> operator
    Options_tree* p = get_single_option(key);
    p->value_iter = p->values.begin();
    return *p;
}

// Specializations for int and double
// If conversion fails, it tryes to extract another type
Options_tree& Options_tree::operator>>(int& val){
    if(value_iter==values.end())
        throw pteros::Pteros_error("No more values!");
    // Set val to current value
    try {
        val = boost::get<int>(*value_iter);
    } catch(boost::bad_get){
        val = boost::get<double>(*value_iter);
    }
    // Advance iterator
    value_iter++;
    return *this;
}

Options_tree& Options_tree::operator>>(double& val){
    if(value_iter==values.end())
        throw pteros::Pteros_error("No more values!");
    // Set val to current value
    try {
        val = boost::get<double>(*value_iter);
    } catch(boost::bad_get){
        val = boost::get<int>(*value_iter);
    }
    // Advance iterator
    value_iter++;
    return *this;
}
