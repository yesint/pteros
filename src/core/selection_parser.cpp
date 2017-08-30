/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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
#include "selection_parser.h"
#include "selection_grammar.h"
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <Eigen/Core>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <unordered_set>
#include <regex>

//-----------------------------------------------------
//  Functions for creating actual selection from AST
//-----------------------------------------------------
using namespace std;
using namespace pteros;
using namespace boost;

#ifdef _DEBUG_PARSER
char* tok_names[] = {
    "TOK_VOID",
    "TOK_MINUS",
    "TOK_UNARY_MINUS",
    "TOK_PLUS",
    "TOK_MULT",
    "TOK_DIV",
    "TOK_POWER",
    "TOK_EQ", // == or =
    "TOK_NEQ", // <> or !=
    "TOK_LT", //<
    "TOK_GT", //>
    "TOK_LEQ", //<=
    "TOK_GEQ", //>=
    // Operations for atom field codes
    "TOK_X",
    "TOK_Y",
    "TOK_Z",
    "TOK_OCC",
    "TOK_BETA",
    // Logic
    "TOK_OR",
    "TOK_AND",
    // Prefixes
    "TOK_NOT",
    "TOK_WITHIN",
    "TOK_SELF",
    //"TOK_PERIODIC",
    //"TOK_OF",
    "TOK_BY",
    "TOK_RESIDUE",
    // text keywords
    "TOK_NAME",
    "TOK_RESNAME",
    "TOK_TAG",
    "TOK_TYPE",
    "TOK_CHAIN",
    // int keywords
    "TOK_RESID",
    "TOK_INDEX",
    "TOK_RESINDEX",
    // all
    "TOK_ALL",
    // Range
    "TOK_TO", // '-' or 'to'
    // Data "TOKens
    "TOK_INT",
    "TOK_UINT",
    "TOK_FLOAT",
    "TOK_STR",
    // Parens
    //"TOK_LPAREN",
    //"TOK_RPAREN",
    // Distances
    //"TOK_DIST",
    "TOK_POINT",
    "TOK_VECTOR",
    "TOK_PLANE",

    "TOK_PRECOMPUTED",
    "TOK_REGEX"
};


void AstNode::dump(int indent){
    for(int i=0;i<indent;++i) cout << '\t';
    cout << tok_names[code] << "{" << endl;
    if(code == TOK_PRECOMPUTED){
        //for(int j=0;j<precomputed.size();++j) cout << precomputed[j]<<" ";
        for(int i=0;i<indent;++i) cout << '\t';
        cout << "\tsz: " << precomputed.size() << endl;
    }

    for(int i=0;i<children.size();++i){        
        AstNode_ptr p;
        try {
            p = boost::get<AstNode_ptr>(children[i]);
            p->dump(indent+1);
        } catch(boost::bad_get) {
            for(int j=0;j<indent;++j) cout << '\t';
            cout << "  " << children[i] << endl;
        }
    }
    for(int i=0;i<indent;++i) cout << '\t';
    cout << "}" << endl;
}

string AstNode::decode(){
    return tok_names[code];
}

#endif

int AstNode::child_as_int(int i){    
    return boost::get<int>(boost::get<AstNode_ptr>(children[i])->children[0]);
}

string AstNode::child_as_str(int i){    
    return boost::get<string>(boost::get<AstNode_ptr>(children[i])->children[0]);
}

float AstNode::child_as_float(int i){        
    return boost::get<float>(boost::get<AstNode_ptr>(children[i])->children[0]);
}

AstNode_ptr& AstNode::child_node(int i){
    return boost::get<AstNode_ptr>(children[i]);
}

bool AstNode::is_coordinate_dependent(){
    if(code == TOK_X || code == TOK_Y || code == TOK_Z || code == TOK_WITHIN
            || code == TOK_POINT || code == TOK_PLANE || code == TOK_VECTOR
      ){
        return true;
    } else {
        return false;
    }
}


Selection_parser::Selection_parser(std::vector<int> *subset):
    has_coord(false),
    starting_subset(subset) { }

Selection_parser::~Selection_parser(){}


bool is_node_pure(AstNode_ptr& node){
    if(node->is_coordinate_dependent()) return false;

    for(int i=0;i<node->children.size();++i){
        try {
            if( is_node_pure(boost::get<AstNode_ptr>(node->children[i])) == false){
                return false;
            }
        } catch(boost::bad_get){};
    }

    return true;
}

void Selection_parser::create_ast(string& sel_str){
#ifdef _DEBUG_PARSER
	cout << "Going to create AST from: " << sel_str <<endl;
#endif
    Grammar g(sel_str);
    tree = g.run();

    if(!is_node_pure(tree)) has_coord = true;
    is_optimized = false; // Not yet optimized

#ifdef _DEBUG_PARSER
    cout << "Is coordinate dependent? " << has_coord << endl;
    tree->dump();
#endif
}



void Selection_parser::do_optimization(AstNode_ptr& node){

    // Skip optimization for trivial terminal nodes
    if(    node->code == TOK_UINT
        || node->code == TOK_INT
        || node->code == TOK_FLOAT
        || node->code == TOK_STR
        || node->code == TOK_REGEX
        || node->code == TOK_X
        || node->code == TOK_Y
        || node->code == TOK_Z
        || node->code == TOK_BETA
        || node->code == TOK_OCC
        || node->code == TOK_TO
        || node->code == TOK_INDEX
        || node->code == TOK_RESINDEX
        || node->code == TOK_RESID
       ) return;

    // Now check if this node does not contain coord-dependent children
    if(is_node_pure(node)){
#ifdef _DEBUG_PARSER
        cout << "Node " << node->decode() << " is pure" << endl;
#endif

        // Node is pure! Check if this is a math expression, which evaluates to constant
        if(    node->code == TOK_PLUS
            || node->code == TOK_MINUS
            || node->code == TOK_MULT
            || node->code == TOK_DIV
            || node->code == TOK_POWER
          )
        {
            // Eval to constant and replace node with float
            float val = get_numeric(node)(0);
            node->code = TOK_FLOAT;
            node->children.clear();
            node->children.push_back(val);
        } else {

            // Not a math expression,so clear all its children and keep precomputed index
            // Set node type to precomputed

            eval_node(node,node->precomputed,nullptr);
            node->children.clear();
            node->code = TOK_PRECOMPUTED;

#ifdef _DEBUG_PARSER
            cout << "Node set to precomputed " << endl;
#endif
        }
    }

    // Optimize AND operations - coord-dependent operand
    // should go second to benefit from subspace optimization
    if(node->code == TOK_AND){
        try{            
            if(  !is_node_pure(node->child_node(0))
               && is_node_pure(node->child_node(1))
               ){

#ifdef _DEBUG_PARSER
                cout << "Node " << node->decode() << " swapped" << endl;
#endif

                node->child_node(0).swap(node->child_node(1));
             }
        } catch (boost::bad_get) {}
    }

    // Go deeper
    for(int i=0;i<node->children.size();++i)
        try {
            do_optimization(node->child_node(i));
        } catch(boost::bad_get){};
}


void Selection_parser::apply(System* system, size_t fr, vector<int>& result){
#ifdef _DEBUG_PARSER
	cout << "Applying to the system with first atom " << system->atoms[0].name << endl;
#endif

    sys = system;
    frame = fr;
    Natoms = sys->num_atoms();

    // For coordinate-dependent selections perform optimization
    // to precompute all pure not-coordinate-depenednt nodes
    if(has_coord && !is_optimized){

#ifdef _DEBUG_PARSER
        cout << "Tree before optimizaton:" << endl;
        tree->dump(0);
#endif
        do_optimization(tree);
        is_optimized = true;

#ifdef _DEBUG_PARSER
        cout << "Tree after optimizaton:" << endl;
        tree->dump(0);
#endif
    }    

    // Eval root node
    eval_node(tree,result,nullptr);

    // Sort result to always get ordered selection index
    //sort(result.begin(),result.end());

#ifdef _DEBUG_PARSER
    int n = 0;
    for(int i=0;i<result.size();++i)
        cout << result[i] << " ";
    cout << endl;    
#endif    
}

void Selection_parser::eval_node(AstNode_ptr& node, vector<int>& result, vector<int>* subspace){
    int i,at,j,k,n;

    // Clear any garbage passed in result
    result.clear();

    //---------------------------------------------------------------------------
    if(node->code == TOK_PRECOMPUTED)
    {
        if(!subspace){
            result = node->precomputed;
        } else {
            // If this node is under subspace we only need to return those
            // atoms, which are in that subspace. Otherwise we can extend
            // the subspace but we don't want this. Thus return intersection:
            std::set_intersection(subspace->begin(),subspace->end(),
                                  node->precomputed.begin(),node->precomputed.end(),
                                  back_inserter(result));
        }        
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_NOT)
    {
        // Logical NOT
        vector<int> res1;
        eval_node(node->child_node(0), res1, subspace);
        // Sort
        //std::sort(res1.begin(),res1.end());
        // res1 is sorted, so we can speed up negation a bit by filling only the "gaps"
        n = res1.size();
        result.reserve(Natoms-n);

        // Special check for empty res1
        if(n==0){
            for(j=0;j<Natoms;++j) result.push_back(j); // All
        } else {

            for(j=0;j<res1[0];++j) result.push_back(j); //Before first
            for(i=1;i<n;++i)
                for(j=res1[i-1]+1;j<res1[i];++j) result.push_back(j); // between any two
            for(j=res1[n-1]+1;j<Natoms;++j) result.push_back(j); // after last
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_OR)
    {
        // Logical OR
        vector<int> res1, res2; // Aux vectors

        eval_node(node->child_node(0), res1, subspace);
        eval_node(node->child_node(1), res2, subspace);

        // Sort
        //std::sort(res1.begin(),res1.end());
        //std::sort(res2.begin(),res2.end());

        std::set_union(res1.begin(),res1.end(),res2.begin(),res2.end(),back_inserter(result));            
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_AND)
    {
        // Logical AND
        vector<int> res1, res2; // Aux vectors

        eval_node(node->child_node(0), res1, subspace);
        // First operand sets a subspace for the second!
        eval_node(node->child_node(1), res2, &res1);

        // Sort
        //std::sort(res1.begin(),res1.end());
        //std::sort(res2.begin(),res2.end());

        std::set_intersection(res1.begin(),res1.end(),res2.begin(),res2.end(),back_inserter(result));
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_NAME)
    {
        int Nchildren = node->children.size(); // Get number of children
        string str;
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            str = node->child_as_str(i);
            if(node->child_node(i)->code == TOK_STR){
                // For normal strings
                if(!subspace){
                    for(at=0;at<Natoms;++at){
                        if(sys->atoms[at].name == str) result.push_back(at);
                    }
                } else {
                    for(j=0;j<subspace->size();++j){
                        at = (*subspace)[j];
                        if(sys->atoms[at].name == str) result.push_back(at);
                    }
                }
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                if(!subspace){
                    for(at=0;at<Natoms;++at){
                        if(std::regex_match(sys->atoms[at].name.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                } else {
                    for(j=0;j<subspace->size();++j){
                        at = (*subspace)[j];
                        if(std::regex_match(sys->atoms[at].name.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                }
            }
        }        
        sort(result.begin(),result.end());
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_RESNAME)
    {
        int Nchildren = node->children.size(); // Get number of children
        string str;
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            str = node->child_as_str(i);
            if(node->child_node(i)->code == TOK_STR){
                // For normal strings
                if(!subspace){
                    for(at=0;at<Natoms;++at)
                        if(sys->atoms[at].resname == str) result.push_back(at);
                } else {
                    for(j=0;j<subspace->size();++j){
                        at = (*subspace)[j];
                        if(sys->atoms[at].resname == str) result.push_back(at);
                    }
                }
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                if(!subspace){
                    for(at=0;at<Natoms;++at){
                        if(std::regex_match(sys->atoms[at].resname.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                } else {
                    for(j=0;j<subspace->size();++j){
                        at = (*subspace)[j];
                        if(std::regex_match(sys->atoms[at].resname.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                }
            }
        }
        sort(result.begin(),result.end());
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_TAG)
    {
        int Nchildren = node->children.size(); // Get number of children
        string str;
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            str = node->child_as_str(i);
            if(node->child_node(i)->code == TOK_STR){
                // For normal strings
                if(!subspace){
                    for(at=0;at<Natoms;++at)
                        if(sys->atoms[at].tag == str) result.push_back(at);
                } else {
                    for(j=0;j<subspace->size();++j){
                        at = (*subspace)[j];
                        if(sys->atoms[at].tag == str) result.push_back(at);
                    }
                }
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                if(!subspace){
                    for(at=0;at<Natoms;++at){
                        if(std::regex_match(sys->atoms[at].tag.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                } else {
                    for(j=0;j<subspace->size();++j){
                        at = (*subspace)[j];
                        if(std::regex_match(sys->atoms[at].tag.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                }
            }
        }
        sort(result.begin(),result.end());
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_CHAIN)
    {
        int Nchildren = node->children.size(); // Get number of children
        char ch;
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            ch = node->child_as_str(i)[0];
            if(!subspace){
                for(at=0;at<Natoms;++at)
                    if(sys->atoms[at].chain == ch) result.push_back(at);
            } else {
                for(j=0;j<subspace->size();++j){
                    at = (*subspace)[j];
                    if(sys->atoms[at].chain == ch) result.push_back(at);
                }
            }
        }
        sort(result.begin(),result.end());
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_RESID)
    {
        int Nchildren = node->children.size(); // Get number of children
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            if(node->child_node(i)->code == TOK_INT) { // Resid could be int!
                k = node->child_as_int(i);
                if(!subspace){
                    for(at=0;at<Natoms;++at)
                        // Even if k is out of range, nothing will crash here
                        if(sys->atoms[at].resid == k) result.push_back(at);
                } else {
                    for(j=0;j<subspace->size();++j){
                        at = (*subspace)[j];
                        // Even if k is out of range, nothing will crash here
                        if(sys->atoms[at].resid == k) result.push_back(at);
                    }
                }
            } else {
                // this is a range, not an integer
                AstNode_ptr range = node->child_node(i);                
                int i1 = range->child_as_int(0);
                int i2 = range->child_as_int(1);
                for(k=i1;k<=i2;++k){
                    if(!subspace){
                        for(at=0;at<Natoms;++at)
                            // Even if k is out of range, nothing will crash here
                            if(sys->atoms[at].resid == k) result.push_back(at);
                    } else {
                        for(j=0;j<subspace->size();++j){
                            at = (*subspace)[j];
                            // Even if k is out of range, nothing will crash here
                            if(sys->atoms[at].resid == k) result.push_back(at);
                        }
                    }
                }
            }
        }
        sort(result.begin(),result.end());
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_RESINDEX)
    {
        int Nchildren = node->children.size(); // Get number of children
        // Cycle over children
        for(i=0;i<Nchildren;++i){                        
            if(node->child_node(i)->code == TOK_UINT) {
                k = node->child_as_int(i);
                if(!subspace){
                    for(at=0;at<Natoms;++at)
                        // Even if k is out of range, nothing will crash here
                        if(sys->atoms[at].resindex == k) result.push_back(at);
                } else {
                    for(j=0;j<subspace->size();++j){
                        at = (*subspace)[j];
                        // Even if k is out of range, nothing will crash here
                        if(sys->atoms[at].resindex == k) result.push_back(at);
                    }
                }
            } else {
                // this is a range, not an integer
                AstNode_ptr range = node->child_node(i);
                int i1 = range->child_as_int(0);
                int i2 = range->child_as_int(1);
                for(k=i1;k<=i2;++k){
                    if(!subspace){
                        for(at=0;at<Natoms;++at)
                            // Even if k is out of range, nothing will crash here
                            if(sys->atoms[at].resindex == k) result.push_back(at);
                    } else {
                        for(j=0;j<subspace->size();++j){
                            at = (*subspace)[j];
                            // Even if k is out of range, nothing will crash here
                            if(sys->atoms[at].resindex == k) result.push_back(at);
                        }
                    }
                }
            }
        }
        sort(result.begin(),result.end());
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_INDEX)
    {
        int Nchildren = node->children.size(); // Get number of children
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            if(node->child_node(i)->code == TOK_UINT) {
                k = node->child_as_int(i);
                // We have to check the range here
                if(k>=0 && k<Natoms)
                    result.push_back(k);
            } else {
                // this is a range, not an integer
                AstNode_ptr range = node->child_node(i);
                int i1 = range->child_as_int(0);
                int i2 = range->child_as_int(1);
                for(k=i1;k<=i2;++k)
                    // We have to check the range here
                    if(k>=0 && k<Natoms)
                        result.push_back(k);
            }
        }
        sort(result.begin(),result.end());
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_WITHIN)
    {
        // Get distance
        double dist = boost::get<float>(node->children[0]);
        // Get PBC
        bool periodic = (boost::get<int>(node->children[2])) ? true : false;
        // Get self
        bool include_self = (boost::get<int>(node->children[3])) ? true : false;               

#ifdef _DEBUG_PARSER
        if(subspace)
            cout << "subspace size: " << subspace->size() << endl;
        else
            cout << "full subspace!" << endl;
#endif
        // Create selections for searching
        Selection dum1(*sys), dum2(*sys);

        // Evaluate enclosed expression
        // Enclosed expression is independent of any subspace!
        // Otherwise the results would be incorrect        
        // Result is returned directly into the index array of selection dum2
        // thus no additional copying
        if(starting_subset && starting_subset->size()){
            eval_node(node->child_node(1), dum2.index, starting_subset);
        } else {
            eval_node(node->child_node(1), dum2.index, nullptr);
        }

        // Prepare selection dum1
        if(!subspace){
            // We are not limited by subspace
            dum1.index.resize(sys->num_atoms());
            for(int i=0;i<sys->num_atoms();++i) dum1.index[i] = i;
        } else {
            // We are limited by subspace
            dum1.index = *subspace;
        }

        // Set frame for both selections
        dum1.set_frame(frame);
        dum2.set_frame(frame);

        search_within(dist,dum1,dum2,result,include_self,periodic);
        // Returned array is sorted already!
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_BY)
    {
        vector<int> res1;
        // Evaluate enclosed expression
        if(starting_subset && starting_subset->size()){
            eval_node(node->child_node(0), res1, starting_subset);
        } else {
            eval_node(node->child_node(0), res1, nullptr);
        }
        int Nsel = res1.size();
        // Select by residue. This respects chain!
        // First make a set of resids we need to search
        std::unordered_set<int> resind;
        for(i=0;i<Nsel;++i){ //over found atoms
            resind.insert(sys->atoms[res1[i]].resindex);
        }

        // Now cycle over all atoms in the system (not a subset!)
        for(at=0;at<Natoms;++at){ // over all atoms
            if(resind.count(sys->atoms[at].resindex)){
                // This resind is needed
                result.push_back(at);
            }
        }
        // Now result is sorted by default here and there are no duplicates
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_ALL)
    {
        result.resize(Natoms);
        for(at=0;at<Natoms;++at) result[at] = at;    
    }

    //---------------------------------------------------------------------------
    // Math logical nodes
    // All of them benefit from subspace optimization
    // All give sorted results

    else if(node->code == TOK_EQ)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!subspace){
            for(at=0;at<Natoms;++at) // over all atoms
                if( op1(at) == op2(at) ) result.push_back(at);
        } else {
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if( op1(at) == op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_NEQ)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!subspace){
            for(at=0;at<Natoms;++at) // over all atoms
                if( op1(at) != op2(at) ) result.push_back(at);
        } else {
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if( op1(at) != op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_LT)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!subspace){
            for(at=0;at<Natoms;++at){ // over all atoms
                if( op1(at) < op2(at) ) result.push_back(at);
            }
        } else {            
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if( op1(at) < op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_GT)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!subspace){
            for(at=0;at<Natoms;++at){ // over all atoms
                if( op1(at) > op2(at) ) result.push_back(at);
            }
        } else {
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if( op1(at) > op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_LEQ)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!subspace){
            for(at=0;at<Natoms;++at){ // over all atoms
                if( op1(at) <= op2(at) ) result.push_back(at);
            }
        } else {
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if( op1(at) <= op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_GEQ)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!subspace){
            for(at=0;at<Natoms;++at){ // over all atoms
                if( op1(at) >= op2(at) ) result.push_back(at);
            }
        } else {
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if( op1(at) >= op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_GEQ)
    {
        throw Pteros_error("Invalid node in the AST!");
    }


}

// Returns callable, which returns value for numeric node for atom at
std::function<float(int)> Selection_parser::get_numeric(AstNode_ptr& node){
    if(node->code == TOK_INT || node->code == TOK_UINT){
        float val = boost::get<int>(node->children[0]);
        return [val](int at){ return val; };
    } else if(node->code == TOK_FLOAT){
        float val = boost::get<float>(node->children[0]);
        return [val](int at){ return val; };
    } else if(node->code == TOK_X){
        return [this](int at){ return sys->traj[frame].coord[at](0); };
    } else if(node->code == TOK_Y){
        return [this](int at){ return sys->traj[frame].coord[at](1); };
    } else if(node->code == TOK_Z){
        return [this](int at){ return sys->traj[frame].coord[at](2); };
    } else if(node->code == TOK_BETA){
        return [this](int at){ return sys->atoms[at].beta; };
    } else if(node->code == TOK_OCC){
        return [this](int at){ return sys->atoms[at].occupancy; };
    } else if(node->code == TOK_INDEX){
        return [](int at){ return at; };
    } else if(node->code == TOK_RESINDEX){
        return [this](int at){ return sys->atoms[at].resindex; };
    } else if(node->code == TOK_RESID){
        return [this](int at){ return sys->atoms[at].resid; };
    } else if(node->code == TOK_UNARY_MINUS){
        auto func = get_numeric(node->child_node(0));
        return [func](int at){ return -func(at); };
    } else if(node->code == TOK_PLUS){        
        auto func1 = get_numeric(node->child_node(0));
        auto func2 = get_numeric(node->child_node(1));
        return [func1,func2](int at){ return func1(at)+func2(at); };
    } else if(node->code == TOK_MINUS){
        auto func1 = get_numeric(node->child_node(0));
        auto func2 = get_numeric(node->child_node(1));
        return [func1,func2](int at){ return func1(at)-func2(at); };
    } else if(node->code == TOK_MULT){
        auto func1 = get_numeric(node->child_node(0));
        auto func2 = get_numeric(node->child_node(1));
        return [func1,func2](int at){ return func1(at)*func2(at); };
    } else if(node->code == TOK_DIV){
        auto func1 = get_numeric(node->child_node(0));
        auto func2 = get_numeric(node->child_node(1));
        return [func1,func2](int at){
            float v = func2(at);
            if(v==0.0) throw Pteros_error("Division by zero in selection!");
            return func1(at)/v;
        };
    } else if(node->code == TOK_POWER) {
        auto func1 = get_numeric(node->child_node(0));
        auto func2 = get_numeric(node->child_node(1));
        return [func1,func2](int at){ return std::pow(func1(at),func2(at)); };
    } else if(node->code == TOK_POINT){
        // Extract point
        Eigen::Vector3f p;

        p(0) = boost::get<float>(node->children[0]);
        p(1) = boost::get<float>(node->children[1]);
        p(2) = boost::get<float>(node->children[2]);

        bool pbc = (boost::get<int>(node->children[3])) ? true : false;

        // Return distance
        if(pbc){
            return [this,p](int at){
                return sys->Box(frame).distance(p, sys->traj[frame].coord[at]);
            };
        } else {
            return [this,&p](int at){
                return (p - sys->traj[frame].coord[at]).norm();
            };
        }

    } else if(node->code == TOK_VECTOR || node->code == TOK_PLANE ){
        // Extract point
        Eigen::Vector3f p;
        p(0) = boost::get<float>(node->children[0]);
        p(1) = boost::get<float>(node->children[1]);
        p(2) = boost::get<float>(node->children[2]);
        // Extract direction vector (or a normal if it's a plane)
        Eigen::Vector3f dir;
        dir(0) = boost::get<float>(node->children[3]);
        dir(1) = boost::get<float>(node->children[4]);
        dir(2) = boost::get<float>(node->children[5]);        

        // pbc
        bool pbc = (boost::get<int>(node->children[6])) ? true : false;

        bool do_plane = (node->code == TOK_PLANE) ? true : false;

        return [this,p,dir,pbc,do_plane](int at){
            Eigen::Vector3f atom = sys->traj[frame].coord[at];

            // Get vector from p to current atom
            Eigen::Vector3f v = atom - p;

            // Project v onto dir
            v = (v.dot(dir)/dir.squaredNorm())*dir;

            if(do_plane){
                // Get closest point on a plane to atom
                v = atom-v;
            } else {
                // Get the end point of projection
                v += p;
            }

            // Return distance between atom and v
            if(pbc){
                return sys->Box(frame).distance(atom, v);
            } else {
                return (atom-v).norm();
            }
        };
    } else {
        throw Pteros_error("Wrong numeric node!");
    }
}
