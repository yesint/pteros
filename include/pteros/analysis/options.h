/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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


#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <vector>

namespace pteros {

class Options;

// Single option such as -t sss fff ggg
class Option {
    friend class Options;
    /// Create options from command line
    /// with nested tasks
    friend void parse_command_line(int argc, char** argv,
                                   Options& toplevel,
                                   std::string task_tag,
                                   std::vector<Options>& tasks);
    /// without nested tasks
    friend void parse_command_line(int argc, char** argv, Options& toplevel);

public:
    /// throw if not a needed type, or more than 1 value
    /// Get single string value
    std::string as_string() const;
    /// Get single int value,
    int as_int() const;
    /// Get single float value
    float as_float() const;
    /// Get single bool value
    bool as_bool() const;

    /// Get list of values of given type
    /// throws if not all elements of needed type
    std::vector<std::string> as_strings() const;
    std::vector<int> as_ints() const;
    std::vector<float> as_floats() const;
    std::vector<bool> as_bools() const;
private:
    // Internal storage is just a vector of strings representing values
    std::vector<std::string> data;
    // names of this option
    std::string name;

    void debug();
};

//---------------------------------------------------------------------------

/// All options obtained from command line including nested options for tasks
class Options {
    /// Create options from command line
    /// with nested tasks
    friend void parse_command_line(int argc, char** argv,
                                   Options& toplevel,
                                   std::string task_tag,
                                   std::vector<Options>& tasks);
    /// without nested tasks
    friend void parse_command_line(int argc, char** argv, Options& toplevel);

public:
    /// Return single option with given name
    const Option& operator()(std::string key) const;
    const Option& operator()(std::string key, std::string default_val);
    bool has(std::string key);
    std::string get_name(){ return task_name; }

    void debug();
private:
    std::vector<Option> data;
    std::string task_name;
};


} // namespace
#endif
