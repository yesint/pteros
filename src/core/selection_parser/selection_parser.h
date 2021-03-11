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
#include <memory>

#include "pteros/core/system.h"
#include "peglib.h"

namespace pteros {

// Custom annoation for peglib ast structure
struct MyAstAnnotation {
    bool is_coord_dependent;
    std::vector<int> precomputed;
    float numeric_value;
};

typedef peg::AstBase<MyAstAnnotation> MyAst;
typedef std::function<void(std::vector<int>&)> result_func_t;

/**
*   Selection parser class.
*   It parses selection text by means of custom recursive-descendent parser
*   The result of parsing in the abstract syntax tree (AST)
*   stored internally in the parser class. The tree is evaluated against the system,
*   which holds the parser.
    This class should never be used directly.
*/
class SelectionParser{
public:
    /** True if there are coordinate keywords in selection.
    *   If true, the parser will persist (not deleted after parsing).
    *   As a result the AST is instantly available. If coordinates change
    *   the AST may be evaluated against the changed coordinates by means
    *   of System.update()
    */
    bool has_coord;  //There are coordinates in selection
    /// Constructor
    SelectionParser(std::vector<int>* subset = nullptr);
    /// Destructor
    virtual ~SelectionParser();
    /// Generates AST from selection string
    void create_ast(std::string& sel_str, System *system);
    /// Apply ast to the given frame. Fills the vector passed from
    /// enclosing System with selection indexes.
    void apply_ast(std::size_t fr, std::vector<int>& result);

private:
    /// AST structure 
    std::shared_ptr<MyAst> tree;

    // AST evaluation stuff
    System* sys;
    int Natoms;
    int frame;    

    void eval_node(const std::shared_ptr<MyAst> &node, std::vector<int>& result);
    std::function<float(int)> get_numeric(const std::shared_ptr<MyAst>& node);
    Eigen::Vector3f get_vector(const std::shared_ptr<MyAst> &node);

    std::vector<int>* starting_subset;
    std::vector<int>* current_subset;

    void precompute(std::shared_ptr<MyAst> &node);
    void optimize_numeric(std::shared_ptr<MyAst> &node);
};

}





