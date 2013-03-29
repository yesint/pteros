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

#include "bindings_options_tree.h"

// For options tree

void Options_tree_from_command_line1(Options_tree* o, boost::python::list& data){
    int n = len(data);
    char* cmd[n];
    for(int i=0;i<n;++i){
        string s = boost::python::extract<string>(data[i]);
        cmd[i] = strdup(s.c_str());
    }
    o->from_command_line(n,cmd);
    // Free allocated memory
    for(int i=0;i<n;++i) free(cmd[i]);
}


int Options_tree_get_value1_int(Options_tree* o, string key){
    return o->get_value<int>(key);
}

int Options_tree_get_value2_int(Options_tree* o, string key, int def){
    return o->get_value<int>(key,def);
}

int Options_tree_get_value1_bool(Options_tree* o, string key){
    return o->get_value<bool>(key);
}

int Options_tree_get_value2_bool(Options_tree* o, string key, int def){
    return o->get_value<bool>(key,def);
}


double Options_tree_get_value1_double(Options_tree* o, string key){
    return o->get_value<double>(key);
}

double Options_tree_get_value2_double(Options_tree* o, string key, double def){
    return o->get_value<double>(key,def);
}

std::string Options_tree_get_value1_string(Options_tree* o, string key){
    return o->get_value<std::string>(key);
}

std::string Options_tree_get_value2_string(Options_tree* o, string key, std::string def){
    return o->get_value<std::string>(key,def);
}


boost::python::list  Options_tree_get_values_int(Options_tree* o, std::string key){
    boost::python::list l;
    std::list<int> r;
    std::list<int>::iterator it;
    r = o->get_values<int>(key);
    for(it=r.begin();it!=r.end();it++) l.append(*it);
    return l;
}

boost::python::list  Options_tree_get_values_double(Options_tree* o, std::string key){
    boost::python::list l;
    std::list<double> r;
    std::list<double>::iterator it;
    r = o->get_values<double>(key);
    for(it=r.begin();it!=r.end();it++) l.append(*it);
    return l;
}

boost::python::list  Options_tree_get_values_string(Options_tree* o, std::string key){
    boost::python::list l;
    std::list<string> r;
    std::list<string>::iterator it;
    r = o->get_values<string>(key);
    for(it=r.begin();it!=r.end();it++) l.append((*it).c_str());
    return l;
}

boost::python::list  Options_tree_get_options(Options_tree* o, std::string key){
    boost::python::list l;
    std::list<Options_tree*> r;
    std::list<Options_tree*>::iterator it;
    r = o->get_options(key);
    for(it=r.begin();it!=r.end();it++) l.append(*it);
    return l;
}


void make_bindings_Options_tree(){

    import_array();

    class_<Options_tree>("Options_tree", init<>())
        .def(init<string>() )
        .def("from_command_line",&Options_tree_from_command_line1)
        .def("to_command_line",&Options_tree::to_command_line)
        .def("to_json_string",&Options_tree::to_json_string)
        .def("from_json_string",&Options_tree::from_json_string)
        .def("get_value_double",&Options_tree_get_value1_double)
        .def("get_value_double",&Options_tree_get_value2_double)
        .def("get_value_int",&Options_tree_get_value1_int)
        .def("get_value_int",&Options_tree_get_value2_int)
        .def("get_value_bool",&Options_tree_get_value1_bool)
        .def("get_value_bool",&Options_tree_get_value2_bool)
        .def("get_value_string",&Options_tree_get_value1_string)
        .def("get_value_string",&Options_tree_get_value2_string)
        .def("get_values_double",&Options_tree_get_values_double)
        .def("get_values_int",&Options_tree_get_values_int)
        .def("get_values_string",&Options_tree_get_values_string)
        .def("count_options",&Options_tree::count_options)
        .def("get_option",&Options_tree::get_option,return_value_policy<reference_existing_object>())
        .def("get_options",&Options_tree_get_options)
    ;
}
