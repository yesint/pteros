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

#include "pteros/core/selection_parser.h"

#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/grid_search.h"
#include <Eigen/Core>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <unordered_set>
#include <regex>

#include "selection_grammar.h"

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


Selection_parser::Selection_parser(){
    has_coord = false;
}

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
            float val = eval_numeric(node,0);
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
    eval_node(tree,result,NULL);

    // Sort result to always get ordered selection index
    sort(result.begin(),result.end());

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

    switch(node->code){
    case  TOK_PRECOMPUTED: {
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
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_NOT: {
        // Logical NOT
        vector<int> res1;
        eval_node(node->child_node(0), res1, subspace);
        // Sort
        std::sort(res1.begin(),res1.end());
        // res1 is sorted, so we can speed up negation a bit by filling only the "gaps"
        n = res1.size();
        for(j=0;j<res1[0];++j) result.push_back(j); //Before first
        for(i=1;i<n;++i)
            for(j=res1[i-1]+1;j<res1[i];++j) result.push_back(j); // between any two
        for(j=res1[n-1]+1;j<Natoms;++j) result.push_back(j); // after last
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_OR: {
        // Logical OR
        vector<int> res1, res2; // Aux vectors

        eval_node(node->child_node(0), res1, subspace);
        eval_node(node->child_node(1), res2, subspace);

        // Sort
        std::sort(res1.begin(),res1.end());
        std::sort(res2.begin(),res2.end());

        std::set_union(res1.begin(),res1.end(),res2.begin(),res2.end(),back_inserter(result));
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_AND: {
        // Logical AND
        vector<int> res1, res2; // Aux vectors

        eval_node(node->child_node(0), res1, subspace);
        // First operand sets a subspace for the second!
        eval_node(node->child_node(1), res2, &res1);

        // Sort
        std::sort(res1.begin(),res1.end());
        std::sort(res2.begin(),res2.end());

        std::set_intersection(res1.begin(),res1.end(),res2.begin(),res2.end(),back_inserter(result));
        return;
    }

    // Coord-independent selections are evaluated only once thus
    // no need to bother with subspace optimizations

    //---------------------------------------------------------------------------
    case  TOK_NAME: {
        int Nchildren = node->children.size(); // Get number of children
        string str;
        // Cycle over children

        for(i=0;i<Nchildren;++i){
            str = node->child_as_str(i);
            if(node->child_node(i)->code == TOK_STR){
                // For normal strings
                for(at=0;at<Natoms;++at){
                    if(sys->atoms[at].name == str) result.push_back(at);
                }
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                for(at=0;at<Natoms;++at){
                    if(std::regex_match(sys->atoms[at].name.c_str(),what,reg)){
                         result.push_back(at);
                    }
                }
            }
        }

        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_RESNAME: {
        int Nchildren = node->children.size(); // Get number of children
        string str;
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            str = node->child_as_str(i);
            if(node->child_node(i)->code == TOK_STR){
                // For normal strings
                for(at=0;at<Natoms;++at)
                    if(sys->atoms[at].resname == str) result.push_back(at);
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                for(at=0;at<Natoms;++at){
                    if(std::regex_match(sys->atoms[at].resname.c_str(),what,reg)){
                         result.push_back(at);
                    }
                }
            }
        }
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_TAG: {
        int Nchildren = node->children.size(); // Get number of children
        string str;
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            str = node->child_as_str(i);
            if(node->child_node(i)->code == TOK_STR){
                // For normal strings
                for(at=0;at<Natoms;++at)
                    if(sys->atoms[at].tag == str) result.push_back(at);
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                for(at=0;at<Natoms;++at){
                    if(std::regex_match(sys->atoms[at].tag.c_str(),what,reg)){
                         result.push_back(at);
                    }
                }
            }
        }
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_CHAIN: {
        int Nchildren = node->children.size(); // Get number of children
        char ch;
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            ch = node->child_as_str(i)[0];
            for(at=0;at<Natoms;++at)
                if(sys->atoms[at].chain == ch) result.push_back(at);
        }
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_RESID: {
        int Nchildren = node->children.size(); // Get number of children
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            if(node->child_node(i)->code == TOK_INT) { // Resid could be int!
                k = node->child_as_int(i);
                for(at=0;at<Natoms;++at)
                    // Even if k is out of range, nothing will crash here
                    if(sys->atoms[at].resid == k) result.push_back(at);
            } else {
                // this is a range, not an integer
                AstNode_ptr range = node->child_node(i);                
                int i1 = range->child_as_int(0);
                int i2 = range->child_as_int(1);
                for(k=i1;k<=i2;++k)
                    for(at=0;at<Natoms;++at)
                        // Even if k is out of range, nothing will crash here
                        if(sys->atoms[at].resid == k) result.push_back(at);
            }
        }
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_RESINDEX: {
        int Nchildren = node->children.size(); // Get number of children
        // Cycle over children
        for(i=0;i<Nchildren;++i){                        
            if(node->child_node(i)->code == TOK_UINT) {
                k = node->child_as_int(i);
                for(at=0;at<Natoms;++at)
                    // Even if k is out of range, nothing will crash here
                    if(sys->atoms[at].resindex == k) result.push_back(at);
            } else {
                // this is a range, not an integer
                AstNode_ptr range = node->child_node(i);
                int i1 = range->child_as_int(0);
                int i2 = range->child_as_int(1);
                for(k=i1;k<=i2;++k)
                    for(at=0;at<Natoms;++at)
                        // Even if k is out of range, nothing will crash here
                        if(sys->atoms[at].resindex == k) result.push_back(at);
            }
        }
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_INDEX: {
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
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_WITHIN: {        
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
        // Enclosed expression is independent on any subspace!
        // Otherwise the results would be incorrect        
        // Result is returned directly into the index array of selection dum2
        // thus no additional copying
        eval_node(node->child_node(1), dum2.index, NULL);

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

        Grid_searcher(dist,dum1,dum2,result,include_self,true,periodic);

        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_BY: {        
        vector<int> res1;
        // Evaluate enclosed expression, ok to pass subspace
        eval_node(node->child_node(0), res1, subspace);
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
        return;
    }
    //---------------------------------------------------------------------------
    case  TOK_ALL:
        result.reserve(Natoms);
        for(at=0;at<Natoms;++at) result.push_back(at);
        return;          
    //---------------------------------------------------------------------------
    // Math logical nodes
    // All of them benefit from subspace optimization
    case  TOK_EQ:
        if(!subspace){
            for(at=0;at<Natoms;++at) // over all atoms
                if(eval_numeric(node->child_node(0),at) ==
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
        } else {
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if(eval_numeric(node->child_node(0),at) ==
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
            }
        }
        return;
    //---------------------------------------------------------------------------
    case  TOK_NEQ:
        if(!subspace){
            for(at=0;at<Natoms;++at) // over all atoms
                if(eval_numeric(node->child_node(0),at) !=
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
        } else {
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if(eval_numeric(node->child_node(0),at) !=
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
            }
        }
        return;
    //---------------------------------------------------------------------------
    case  TOK_LT:
        if(!subspace){            
            for(at=0;at<Natoms;++at){ // over all atoms
                if(eval_numeric(node->child_node(0),at) <
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
            }
        } else {            
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if(eval_numeric(node->child_node(0),at) <
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
            }
        }
        return;
    //---------------------------------------------------------------------------
    case  TOK_GT:        
        if(!subspace){            
            for(at=0;at<Natoms;++at) // over all atoms
                if(eval_numeric(node->child_node(0),at) >
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);            
        } else {            
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if(eval_numeric(node->child_node(0),at) >
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
            }
        }
        return;
    //---------------------------------------------------------------------------
    case  TOK_LEQ:
        if(!subspace){
            for(at=0;at<Natoms;++at) // over all atoms
                if(eval_numeric(node->child_node(0),at) <=
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
        } else {
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if(eval_numeric(node->child_node(0),at) <=
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
            }
        }
        return;
    //---------------------------------------------------------------------------
    case  TOK_GEQ:
        if(!subspace){
            for(at=0;at<Natoms;++at) // over all atoms
                if(eval_numeric(node->child_node(0),at) >=
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
        } else {
            for(int i=0;i<subspace->size();++i){ // over subspace
                at = (*subspace)[i];
                if(eval_numeric(node->child_node(0),at) >=
                   eval_numeric(node->child_node(1),at)
                  ) result.push_back(at);
            }
        }
        return;

    }   
}

float Selection_parser::eval_numeric(AstNode_ptr& node, int at){    
    if(node->code == TOK_INT){
        return boost::get<int>(node->children[0]);
    } else if(node->code == TOK_UINT){
        return boost::get<int>(node->children[0]);
    } else if(node->code == TOK_FLOAT){        
        return boost::get<float>(node->children[0]);
    } else if(node->code == TOK_X){        
        return sys->traj[frame].coord[at](0);
    } else if(node->code == TOK_Y){
        return sys->traj[frame].coord[at](1);
    } else if(node->code == TOK_Z){        
        return sys->traj[frame].coord[at](2);
    } else if(node->code == TOK_BETA){
        return sys->atoms[at].beta;
    } else if(node->code == TOK_OCC){
        return sys->atoms[at].occupancy;
    } else if(node->code == TOK_INDEX){
        return at;
    } else if(node->code == TOK_RESINDEX){
        return sys->atoms[at].resindex;
    } else if(node->code == TOK_UNARY_MINUS){
        return -eval_numeric(node->child_node(0),at);
    } else if(node->code == TOK_PLUS){
        return eval_numeric(node->child_node(0),at)
             + eval_numeric(node->child_node(1),at);
    } else if(node->code == TOK_MINUS){
        return eval_numeric(node->child_node(0),at)
             - eval_numeric(node->child_node(1),at);
    } else if(node->code == TOK_MULT){
        return eval_numeric(node->child_node(0),at)
             * eval_numeric(node->child_node(1),at);
    } else if(node->code == TOK_DIV){
        float v = eval_numeric(node->child_node(0),at);
        if(v==0.0) throw Pteros_error("Divition by zero in selection!");
        return eval_numeric(node->child_node(1),at) / v;
    } else if(node->code == TOK_POWER) {
        return std::pow( eval_numeric(node->child_node(0),at),
                         eval_numeric(node->child_node(1),at) );
    } else if(node->code == TOK_POINT){
        // Extract point
        Eigen::Vector3f p;

        p(0) = eval_numeric(node->child_node(0),at);
        p(1) = eval_numeric(node->child_node(1),at);
        p(2) = eval_numeric(node->child_node(2),at);

        bool pbc = (boost::get<int>(node->children[3])) ? true : false;

        // Return distance
        if(pbc){
            return sys->Box(frame).distance(p, sys->traj[frame].coord[at]);
        } else {
            return (p - sys->traj[frame].coord[at]).norm();
        }

    } else if(node->code == TOK_VECTOR || node->code == TOK_PLANE ){
        // Extract point
        Eigen::Vector3f p;
        p(0) = eval_numeric(node->child_node(0),at);
        p(1) = eval_numeric(node->child_node(1),at);
        p(2) = eval_numeric(node->child_node(2),at);
        // Extract direction vector (or a normal if it's a plane)
        Eigen::Vector3f dir;
        dir(0) = eval_numeric(node->child_node(3),at);
        dir(1) = eval_numeric(node->child_node(4),at);
        dir(2) = eval_numeric(node->child_node(5),at);

        // pbc
        bool pbc = (boost::get<int>(node->children[6])) ? true : false;

        Eigen::Vector3f atom = sys->traj[frame].coord[at];

        // Get vector from p to current atom
        Eigen::Vector3f v = atom - p;

        // Project v onto dir
        v = (v.dot(dir)/dir.squaredNorm())*dir;

        if(node->code == TOK_PLANE){
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
    }
}
