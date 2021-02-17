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

#pragma once

#include <string>
#include <vector>
#include <Eigen/Core>

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
    /// Throws if not a needed type, or more than 1 value
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
    Eigen::VectorXf as_VectorXf() const;
private:
    // Internal storage is just a vector of strings representing values
    std::vector<std::string> data;
    // name of this option
    std::string name;
#ifdef DEBUG
    void debug();
#endif
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
#ifdef DEBUG
    void debug();
#endif
private:
    std::vector<Option> data;
    std::string task_name;
};


} // namespace




