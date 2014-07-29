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

#define _DEBUG_PARSER

#ifndef SELECTION_PARSER_H
#define SELECTION_PARSER_H

#include <string>
#include <vector>
#include <memory>

#include "pteros/core/system.h"
#include <boost/variant.hpp>

namespace pteros {

// Codes of tokens
enum Codes {
    //TOK_VOID,
    TOK_MINUS,
    TOK_UNARY_MINUS,
    TOK_PLUS,
    TOK_MULT,
    TOK_DIV,
    TOK_POWER,
    TOK_EQ, // == or =
    TOK_NEQ, // <> or !=
    TOK_LT, //<
    TOK_GT, //>
    TOK_LEQ, //<=
    TOK_GEQ, //>=
    // Operations for atom field codes
    TOK_X,
    TOK_Y,
    TOK_Z,
    TOK_OCC,
    TOK_BETA,
    // Logic
    TOK_OR,
    TOK_AND,
    // Prefixes
    TOK_NOT,
    TOK_WITHIN,    
    //TOK_PERIODIC,
    //TOK_OF,
    TOK_BY,
    TOK_RESIDUE,
    // text keywords
    TOK_NAME,
    TOK_RESNAME,
    TOK_TAG,
    TOK_TYPE,
    TOK_CHAIN,
    // int keywords
    TOK_RESID,
    TOK_INDEX,
    TOK_RESINDEX,
    // all
    TOK_ALL,
    // Range
    TOK_TO, // '-' or 'to'
    // Data tokens
    TOK_INT,
    TOK_UINT,
    TOK_FLOAT,
    TOK_STR,
    // Parens
    //TOK_LPAREN,
    //TOK_RPAREN,
    // Distances
    //TOK_DIST,
    TOK_POINT,
    TOK_VECTOR,
    TOK_PLANE,

    TOK_PRECOMPUTED,
    TOK_REGEX
};

struct AstNode; // Forward declaration
// An element of the tree is either a recursive sub-tree or a leaf
typedef boost::variant<
    float,
    char,
    int,   
    std::string,
    std::shared_ptr<AstNode>
> ast_element;

// The tree itself
struct AstNode {
    Codes code; //Code of operation
    std::vector<ast_element> children;
    std::vector<int> precomputed; // Precomputed indexes for coordinate-independent nodes

    bool is_coordinate_dependent();
    // Returns child elements
    int child_as_int(int i);
    char child_as_char(int i);
    float child_as_float(int i);
    bool child_as_bool(int i);
    float child_as_float_or_int(int i);
    std::string child_as_str(int i);
    std::shared_ptr<AstNode>& child_node(int i);


#ifdef _DEBUG_PARSER
    void dump(int indent=0);
    std::string decode();
#endif
};

typedef std::shared_ptr<AstNode> AstNode_ptr;

/**
*   Selection parser class.
*   It parses selection text by means of custom recursive-descendent parser
*   The result of parsing in the abstract syntax tree (AST)
*   stored internally in the parser class. The tree is evaluated against the system,
*   which holds the parser.
    This class should never be used directly.
*/

class Selection_parser{    
public:
    /** True if there are coordinate keywords in selection.
    *   If true, the parser will persist (not deleted after parsing).
    *   As a result the AST is instantly available. If coordinates change
    *   the AST may be evaluated against the changed coordinates by means
    *   of System.update()
    */
    bool has_coord;  //There are coordinates in selection
    /// Constructor
    Selection_parser();
    /// Destructor
    virtual ~Selection_parser();
    /// Generates AST from selection string
    void create_ast(std::string&);

    /// Apply ast to the system. Fills the vector passed from
    /// enclosing System with selection indexes.
    void apply(System* system, std::size_t fr, std::vector<int>& result);

private:
    /// AST structure
    std::shared_ptr<AstNode> tree;

    // AST evaluation stuff
    System* sys;
    int Natoms;
    int frame;

    void eval_node(AstNode_ptr& node, std::vector<int>& result, std::vector<int>* subspace);
    float eval_numeric(AstNode_ptr& node, int at);
    void do_optimization(AstNode_ptr& node);

    bool is_optimized;        
};

}
#endif /* SELECTION_PARSER_H */
