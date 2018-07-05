/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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

#include "selection_parser.h"

#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <Eigen/Core>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <unordered_set>
#include <regex>

using namespace std;
using namespace pteros;
using namespace boost;


//===============================================

// We derive from normal peg::parser
class Pteros_PEG_parser: public peg::parser {
public:
    Pteros_PEG_parser(const char* s): peg::parser(s) {
        enable_ast();
        enable_packrat_parsing();
        log = [&](size_t ln, size_t col, const string& msg) {
            error_message = fmt::format("{}:{}",col,msg);
        };
    }

    virtual ~Pteros_PEG_parser(){}

    string error_message;
};


// Instance of the parser itself
Pteros_PEG_parser _parser(R"(
        LOGICAL_EXPR       <-  LOGICAL_OPERAND (LOGICAL_OPERATOR LOGICAL_OPERAND)*
        LOGICAL_OPERATOR   <-  'or' / 'and'
        LOGICAL_OPERAND    <-  (NOT / BYRES)? ( '(' LOGICAL_EXPR ')' / ALL / NUM_COMPARISON / KEYWORD_EXPR / WITHIN )
        ALL                <-  'all'
        NOT                <-  'not'
        BYRES              <-  'by' 'residue' / 'same' 'residue' 'as'

        NUM_COMPARISON     <-  NUM_EXPR COMPARISON_OPERATOR NUM_EXPR (COMPARISON_OPERATOR NUM_EXPR)?
        COMPARISON_OPERATOR <- < '<' / '>' / '=' / '==' / '<=' / '>=' / '<>' / '!=' >
        NUM_EXPR           <- NUM_TERM (PLUS_MINUS NUM_TERM)*
        NUM_TERM           <- NUM_POWER (DIV_MUL NUM_POWER)*
        NUM_POWER          <- NUM_FACTOR (POW NUM_FACTOR)?
        NUM_FACTOR         <- UNARY_MINUS? ( '(' NUM_EXPR ')' / X / Y / Z / BETA / OCC / RESINDEX / INDEX / RESID / DIST) / FLOAT
        PLUS_MINUS         <- < '+' / '-' >
        DIV_MUL            <- < '*' / '/' >
        POW                <- < '^' / '**' >
        UNARY_MINUS        <- < '-' >

        X                  <- < 'X' / 'x' >
        Y                  <- < 'Y' / 'y' >
        Z                  <- < 'Z' / 'z' >
        BETA               <- < 'beta' >
        OCC                <- < 'occupancy' / 'occ' >
        RESINDEX           <- < 'resindex' >
        INDEX              <- < 'index' >
        RESID              <- < 'resid' >

        DIST               <- ('dist' / 'distance') (POINT / VECTOR / PLANE)
        POINT              <- 'point' PBC? FLOAT FLOAT FLOAT
        VECTOR             <- 'vector' PBC? FLOAT FLOAT FLOAT FLOAT FLOAT FLOAT
        PLANE              <- 'plane' PBC? FLOAT FLOAT FLOAT FLOAT FLOAT FLOAT

        FLOAT              <- < INTEGER ('.' [0-9]+ )? ( ('e' / 'E' ) INTEGER )? >
        INTEGER            <- < ('-' / '+')? [0-9]+ >

        KEYWORD_EXPR       <- STR_KEYWORD_EXPR / INT_KEYWORD_EXPR
        STR_KEYWORD_EXPR   <- STR_KEYWORD (STR / REGEX)+
        STR_KEYWORD        <- < 'name' / 'resname' / 'tag' / 'chain' / 'type' >
        STR                <- !('or'/'and') < [a-zA-Z0-9]+ >

        INT_KEYWORD_EXPR   <- INT_KEYWORD (RANGE / INTEGER)+
        INT_KEYWORD        <- < 'index' / 'resindex' / 'resid' >
        RANGE              <- INTEGER ('-'/'to'/':') INTEGER

        WITHIN             <- 'within' FLOAT ((PBC SELF)? / (SELF PBC)? / PBC? / SELF?) 'of' LOGICAL_OPERAND
        PBC                <- < 'pbc' / 'nopbc' / 'periodic' / 'nonperiodic' >
        SELF               <- < 'self' / 'noself' >

        %whitespace      <-  [ \t\r\n]*

        REGEX              <- '"' <(!'"' .)*> '"' / "'" <(!"'" .)*> "'"
    )");



//===============================================


bool is_node_coordinate_dependent(const std::shared_ptr<peg::Ast>& node){
    if(node->name == "X" || node->name == "Y" || node->name == "Z" || node->name == "WITHIN"
            || node->name == "POINT" || node->name == "PLANE" || node->name == "VECTOR"
      ){
        return true;
    } else {
        return false;
    }
}


Selection_parser::Selection_parser(std::vector<int> *subset):
    has_coord(false),
    starting_subset(subset)
{

}

Selection_parser::~Selection_parser(){}


bool is_node_pure(const std::shared_ptr<peg::Ast>& node){
    if(is_node_coordinate_dependent(node)) return false;

    for(int i=0;i<node->nodes.size();++i){
        if(is_node_pure(node->nodes[i]) == false) return false;
    }

    return true;
}

void Selection_parser::create_ast(string& sel_str){            
    if (_parser.parse(sel_str.c_str(), tree)) {
        tree = peg::AstOptimizer(true).optimize(tree);

        cout << peg::ast_to_s(tree);

    } else {
        throw Pteros_error(_parser.error_message);
    }

    if(!is_node_pure(tree)) has_coord = true;
    is_optimized = false; // Not yet optimized
}



void Selection_parser::do_optimization(std::shared_ptr<peg::Ast> &node){        

    // Skip optimization for trivial terminal nodes
    if(    node->name == "INTEGER"
        || node->name == "FLOAT"
        || node->name == "STR"
        || node->name == "REGEX"
        || node->name == "X"
        || node->name == "Y"
        || node->name == "Z"
        || node->name == "BETA"
        || node->name == "OCC"
        || node->name == "RANGE"
        || node->name == "INDEX"
        || node->name == "RESINDEX"
        || node->name == "RESID"
       ) return;
/*
    // Now check if this node does not contain coord-dependent children
    if(is_node_pure(node)){
        // Node is pure! Check if this is a math expression, which evaluates to constant
        if(    node->name == TOK_PLUS
            || node->name == TOK_MINUS
            || node->name == TOK_MULT
            || node->name == TOK_DIV
            || node->name == TOK_POWER
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
        }
    }

    // Optimize AND operations - coord-dependent operand
    // should go second to benefit from subspace optimization
    if(node->code == TOK_AND){
        try{            
            if(  !is_node_pure(node->child_node(0))
               && is_node_pure(node->child_node(1))
               ){

                node->child_node(0).swap(node->child_node(1));
             }
        } catch (boost::bad_get) {}
    }
    */

    // Go deeper
    for(int i=0;i<node->nodes.size();++i) do_optimization(node->nodes[i]);

}


void Selection_parser::apply(System* system, size_t fr, vector<int>& result){
    sys = system;
    frame = fr;
    Natoms = sys->num_atoms();

    // For coordinate-dependent selections perform optimization
    // to precompute all pure not-coordinate-depenednt nodes
    if(has_coord && !is_optimized){
        do_optimization(tree);
        is_optimized = true;
    }    

    // Eval root node
    eval_node(tree,result,nullptr);
}

void Selection_parser::eval_node(const std::shared_ptr<peg::Ast> &node, vector<int>& result, vector<int>* subspace){

    int i,at,j,k,n;

    // Clear any garbage passed in result
    result.clear();

    // Pointer to any restricting set (either subspace or starting subset)
    vector<int>* restr;
    if(subspace){
        restr = subspace;
    } else if(starting_subset) {
        restr = starting_subset;
    } else {
        restr = nullptr;
    }        

    if(node->name == "NUM_COMPARISON")
    {
        auto op  = node->nodes[1]->token;
        bool pure1,pure2;
        auto op1 = get_numeric(node->nodes[0],pure1);
        auto op2 = get_numeric(node->nodes[2],pure2);

        // Function to evaluate
        std::function<bool(float,float)> comparison;

        if(op == "=" || op== "=="){
            comparison = [](float a, float b){ return a==b; };
        } else if (op == "!=" || op=="<>"){
            comparison = [](float a, float b){ return a!=b; };
        } else if (op == "<"){
            comparison = [](float a, float b){ return a<b; };
        } else if (op == ">"){
            comparison = [](float a, float b){ return a>b; };
        } else if (op == "<="){
            comparison = [](float a, float b){ return a<=b; };
        } else if (op == ">="){
            comparison = [](float a, float b){ return a>=b; };
        }

        // If both operands are pure evaluate function in place
        // for atom 0 (it does not matter which one to use)
        if(pure1 && pure2){
            LOG()->warn("Meaningless expression in selection");
            if(!comparison(op1(0),op2(0))) throw Pteros_error("False arithmetic comparison");
        } else {
            if(!restr){
                for(at=0;at<Natoms;++at) // over all atoms
                    if( comparison(op1(at),op2(at)) ) result.push_back(at);
            } else {
                for(int i=0;i<restr->size();++i){ // over restr
                    at = (*restr)[i];
                    if( comparison(op1(at),op2(at)) ) result.push_back(at);
                }
            }
        }
    }
/*
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
        eval_node(node->child_node(0), res1, restr);

        n = res1.size();        

        // Special check for empty res1
        if(n==0){
            for(j=0;j<Natoms;++j) result.push_back(j); // All
        } else if(!restr) {
            // Without subset
            result.reserve(Natoms-n);
            for(j=0;j<res1[0];++j) result.push_back(j); //Before first
            for(i=1;i<n;++i)
                for(j=res1[i-1]+1;j<res1[i];++j) result.push_back(j); // between any two
            for(j=res1[n-1]+1;j<Natoms;++j) result.push_back(j); // after last
        } else {
            // With restriction
            std::set_difference(restr->begin(),restr->end(), res1.begin(),res1.end(), back_inserter(result));
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_OR)
    {
        // Logical OR
        vector<int> res1, res2; // Aux vectors

        eval_node(node->child_node(0), res1, restr);
        eval_node(node->child_node(1), res2, restr);

        std::set_union(res1.begin(),res1.end(),res2.begin(),res2.end(),back_inserter(result));            
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_AND)
    {
        // Logical AND
        vector<int> res1, res2; // Aux vectors

        eval_node(node->child_node(0), res1, restr);
        // First operand sets a restr for the second!
        eval_node(node->child_node(1), res2, &res1);

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
                if(!restr){
                    for(at=0;at<Natoms;++at){
                        if(sys->atoms[at].name == str) result.push_back(at);
                    }
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
                        if(sys->atoms[at].name == str) result.push_back(at);
                    }
                }
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                if(!restr){
                    for(at=0;at<Natoms;++at){
                        if(std::regex_match(sys->atoms[at].name.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
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
    else if(node->code == TOK_TYPE)
    {
        if(!sys->force_field.ready) throw Pteros_error("Can't select type names from system with no valid force field!");

        int Nchildren = node->children.size(); // Get number of children
        string str;
        // Cycle over children
        for(i=0;i<Nchildren;++i){
            str = node->child_as_str(i);
            if(node->child_node(i)->code == TOK_STR){
                // For normal strings
                if(!restr){
                    for(at=0;at<Natoms;++at){
                        if(sys->atoms[at].type_name == str) result.push_back(at);
                    }
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
                        if(sys->atoms[at].type_name == str) result.push_back(at);
                    }
                }
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                if(!restr){
                    for(at=0;at<Natoms;++at){
                        if(std::regex_match(sys->atoms[at].type_name.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
                        if(std::regex_match(sys->atoms[at].type_name.c_str(),what,reg)){
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
                if(!restr){
                    for(at=0;at<Natoms;++at)
                        if(sys->atoms[at].resname == str) result.push_back(at);
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
                        if(sys->atoms[at].resname == str) result.push_back(at);
                    }
                }
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                if(!restr){
                    for(at=0;at<Natoms;++at){
                        if(std::regex_match(sys->atoms[at].resname.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
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
                if(!restr){
                    for(at=0;at<Natoms;++at)
                        if(sys->atoms[at].tag == str) result.push_back(at);
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
                        if(sys->atoms[at].tag == str) result.push_back(at);
                    }
                }
            } else if(node->child_node(i)->code == TOK_REGEX){
                // For regex
                std::cmatch what;
                std::regex reg(str);
                if(!restr){
                    for(at=0;at<Natoms;++at){
                        if(std::regex_match(sys->atoms[at].tag.c_str(),what,reg)){
                             result.push_back(at);
                        }
                    }
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
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
            if(!restr){
                for(at=0;at<Natoms;++at)
                    if(sys->atoms[at].chain == ch) result.push_back(at);
            } else {
                for(j=0;j<restr->size();++j){
                    at = (*restr)[j];
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
                if(!restr){
                    for(at=0;at<Natoms;++at)
                        // Even if k is out of range, nothing will crash here
                        if(sys->atoms[at].resid == k) result.push_back(at);
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
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
                    if(!restr){
                        for(at=0;at<Natoms;++at)
                            // Even if k is out of range, nothing will crash here
                            if(sys->atoms[at].resid == k) result.push_back(at);
                    } else {
                        for(j=0;j<restr->size();++j){
                            at = (*restr)[j];
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
                if(!restr){
                    for(at=0;at<Natoms;++at)
                        // Even if k is out of range, nothing will crash here
                        if(sys->atoms[at].resindex == k) result.push_back(at);
                } else {
                    for(j=0;j<restr->size();++j){
                        at = (*restr)[j];
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
                    if(!restr){
                        for(at=0;at<Natoms;++at)
                            // Even if k is out of range, nothing will crash here
                            if(sys->atoms[at].resindex == k) result.push_back(at);
                    } else {
                        for(j=0;j<restr->size();++j){
                            at = (*restr)[j];
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
            eval_node(node->child_node(1), dum2._index, starting_subset);
        } else {
            eval_node(node->child_node(1), dum2._index, nullptr);
        }

        // Prepare selection dum1
        if(!subspace){
            // We are not limited by subspace
            dum1._index.resize(sys->num_atoms());
            for(int i=0;i<sys->num_atoms();++i) dum1._index[i] = i;
        } else {
            // We are limited by subspace
            dum1._index = *subspace;
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
        if(!restr){
            result.resize(Natoms);
            for(at=0;at<Natoms;++at) result[at] = at;
        } else {
            result = *restr; // subspace optimization
        }
    }

    //---------------------------------------------------------------------------
    // Math logical nodes
    // All of them benefit from restr optimization
    // All give sorted results

    else if(node->code == TOK_EQ)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!restr){
            for(at=0;at<Natoms;++at) // over all atoms
                if( op1(at) == op2(at) ) result.push_back(at);
        } else {
            for(int i=0;i<restr->size();++i){ // over restr
                at = (*restr)[i];
                if( op1(at) == op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_NEQ)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!restr){
            for(at=0;at<Natoms;++at) // over all atoms
                if( op1(at) != op2(at) ) result.push_back(at);
        } else {
            for(int i=0;i<restr->size();++i){ // over restr
                at = (*restr)[i];
                if( op1(at) != op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_LT)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!restr){
            for(at=0;at<Natoms;++at){ // over all atoms
                if( op1(at) < op2(at) ) result.push_back(at);
            }
        } else {            
            for(int i=0;i<restr->size();++i){ // over restr
                at = (*restr)[i];
                if( op1(at) < op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_GT)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!restr){
            for(at=0;at<Natoms;++at){ // over all atoms
                if( op1(at) > op2(at) ) result.push_back(at);
            }
        } else {
            for(int i=0;i<restr->size();++i){ // over restr
                at = (*restr)[i];
                if( op1(at) > op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_LEQ)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!restr){
            for(at=0;at<Natoms;++at){ // over all atoms
                if( op1(at) <= op2(at) ) result.push_back(at);
            }
        } else {
            for(int i=0;i<restr->size();++i){ // over restr
                at = (*restr)[i];
                if( op1(at) <= op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_GEQ)
    {
        auto op1 = get_numeric(node->child_node(0));
        auto op2 = get_numeric(node->child_node(1));
        if(!restr){
            for(at=0;at<Natoms;++at){ // over all atoms
                if( op1(at) >= op2(at) ) result.push_back(at);
            }
        } else {
            for(int i=0;i<restr->size();++i){ // over restr
                at = (*restr)[i];
                if( op1(at) >= op2(at) ) result.push_back(at);
            }
        }
    }

    //---------------------------------------------------------------------------
    else if(node->code == TOK_GEQ)
    {
        throw Pteros_error("Invalid node in the AST!");
    }

*/
}

// Returns callable, which returns value for numeric node for atom at
std::function<float(int)> Selection_parser::get_numeric(const std::shared_ptr<peg::Ast> &node, bool& is_pure){

    std::function<float(int)> res;
    is_pure = true;

    // terminals
    if(node->name == "INTEGER"){
        float val = stol(node->token);
        res = [val](int at){ return val; };
    } else if(node->name == "FLOAT"){
        float val = stof(node->token);
        res = [val](int at){ return val; };
    } else if(node->name == "X"){
        res =[this](int at){ return sys->traj[frame].coord[at](0); };
        is_pure = false;
    } else if(node->name == "Y"){
        res = [this](int at){ return sys->traj[frame].coord[at](1); };
        is_pure = false;
    } else if(node->name == "Z"){
        res = [this](int at){ return sys->traj[frame].coord[at](2); };
        is_pure = false;
    } else if(node->name == "BETA"){
        res = [this](int at){ return sys->atoms[at].beta; };
    } else if(node->name == "OCC"){
        res = [this](int at){ return sys->atoms[at].occupancy; };
    } else if(node->name == "INDEX"){
        res = [](int at){ return at; };
    } else if(node->name == "RESINDEX"){
        res = [this](int at){ return sys->atoms[at].resindex; };
    } else if(node->name == "RESID"){
        res = [this](int at){ return sys->atoms[at].resid; };

    // Compounds
    } else if(node->name == "UNARY_MINUS"){
        auto func = get_numeric(node->nodes[0],is_pure);
        res = [func](int at){ return -func(at); };
    } else if(node->name == "NUM_EXPR" || node->name == "NUM_TERM"){
        bool p1,p2;
        auto func1 = get_numeric(node->nodes[0],p1);
        auto op = node->nodes[1]->token;
        auto func2 = get_numeric(node->nodes[2],p2);
        is_pure = p1 && p2;

        if     (op=="+")
            res = [func1,func2](int at){ return func1(at)+func2(at); };
        else if(op=="-")
            res = [func1,func2](int at){ return func1(at)-func2(at); };
        else if(op=="*")
            res = [func1,func2](int at){ return func1(at)*func2(at); };
        else if(op=="/") {
            res = [func1,func2](int at){
                float v = func2(at);
                if(v==0.0) throw Pteros_error("Division by zero in selection!");
                return func1(at)/v;
            };
        }
    }

    return res;

    /*
    else if(node->name == TOK_MINUS){
        auto func1 = get_numeric(node->child_node(0));
        auto func2 = get_numeric(node->child_node(1));
        return [func1,func2](int at){ return func1(at)-func2(at); };
    } else if(node->name == TOK_MULT){
        auto func1 = get_numeric(node->child_node(0));
        auto func2 = get_numeric(node->child_node(1));
        return [func1,func2](int at){ return func1(at)*func2(at); };
    } else if(node->name == TOK_DIV){
        auto func1 = get_numeric(node->child_node(0));
        auto func2 = get_numeric(node->child_node(1));
        return [func1,func2](int at){
            float v = func2(at);
            if(v==0.0) throw Pteros_error("Division by zero in selection!");
            return func1(at)/v;
        };
    } else if(node->name == TOK_POWER) {
        auto func1 = get_numeric(node->child_node(0));
        auto func2 = get_numeric(node->child_node(1));
        return [func1,func2](int at){ return std::pow(func1(at),func2(at)); };
    } else if(node->name == TOK_POINT){
        // Extract point
        Eigen::Vector3f p;

        p(0) = boost::get<float>(node->children[0]);
        p(1) = boost::get<float>(node->children[1]);
        p(2) = boost::get<float>(node->children[2]);

        bool pbc = (boost::get<int>(node->children[3])) ? true : false;

        // Return distance
        if(pbc){
            return [this,p](int at){
                return sys->box(frame).distance(p, sys->traj[frame].coord[at]);
            };
        } else {
            return [this,&p](int at){
                return (p - sys->traj[frame].coord[at]).norm();
            };
        }

    } else if(node->name == TOK_VECTOR || node->name == TOK_PLANE ){
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

        bool do_plane = (node->name == TOK_PLANE) ? true : false;

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
                return sys->box(frame).distance(atom, v);
            } else {
                return (atom-v).norm();
            }
        };
    } else {
        throw Pteros_error("Wrong numeric node!");        
    }
    */
}


