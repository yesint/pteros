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
#include <boost/range/counting_range.hpp>
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
        enable_ast<MyAst>();
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
        LOGICAL_OPERATOR   <-  < 'or' / 'and' >
        LOGICAL_OPERAND    <-  (NOT / BYRES)? ( '(' LOGICAL_EXPR ')' / ALL / NUM_COMPARISON / KEYWORD_EXPR / WITHIN )
        ALL                <-  < 'all' >
        NOT                <-  < 'not' >
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

        WITHIN             <- 'within' FLOAT (PBC SELF / SELF PBC / PBC / SELF)? 'of' LOGICAL_OPERAND

        PBC                <- < 'pbc' / 'nopbc' / 'periodic' / 'nonperiodic' >
        SELF               <- < 'self' / 'noself' >

        %whitespace      <-  [ \t\r\n]*

        REGEX              <- '"' <(!'"' .)*> '"' / "'" <(!"'" .)*> "'"
    )");



//===============================================


bool is_node_coordinate_dependent(const std::shared_ptr<MyAst>& node){
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
    starting_subset(subset),
    sys(nullptr)
{

}

Selection_parser::~Selection_parser(){}


void set_coord_dependence(const std::shared_ptr<MyAst>& node){
    node->is_coord_dependent = is_node_coordinate_dependent(node);

    for(int i=0;i<node->nodes.size();++i){
        set_coord_dependence(node->nodes[i]);
    }
}

void Selection_parser::create_ast(string& sel_str){            
    if (_parser.parse(sel_str.c_str(), tree)) {
        tree = peg::AstOptimizer(true).optimize(tree);

        //cout << peg::ast_to_s(tree);

        set_coord_dependence(tree);
    } else {
        throw Pteros_error(_parser.error_message);
    }

    if(tree->is_coord_dependent) has_coord = true;
}

void Selection_parser::apply(System* system, size_t fr, vector<int>& result){    
    frame = fr;    

    if(starting_subset)
        current_subset = starting_subset;
    else
        current_subset = nullptr;

    // If this is first usage ever or usage for new system
    // make new evaluation function
    if(!sys || sys!=system){
        // Set new system
        sys = system;
        Natoms = sys->num_atoms();
        // generate new evaluation function
        evaluation_function = eval_node(tree);
    }

    // Apply evaluation function
    evaluation_function(result);
}

result_func_t Selection_parser::eval_node(const std::shared_ptr<MyAst> &node){

    result_func_t returned_function;

    // Here starts evaluation

    if(node->name == "NUM_COMPARISON")
    {
        auto op  = node->nodes[1]->token;
        bool pure1,pure2;
        auto op1 = get_numeric(node->nodes[0]);
        auto op2 = get_numeric(node->nodes[2]);

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

        if(!node->nodes[0]->is_coord_dependent && !node->nodes[2]->is_coord_dependent){
            // If both operands are pure return empty vector
            LOG()->warn("Meaningless expression in selection");
            if(!comparison(op1(0),op2(0))) throw Pteros_error("False arithmetic comparison");

            // Return empty
            returned_function = [](vector<int>& result){
                result.clear();
            };

        } else {            

            returned_function = [this,op1,op2,comparison](vector<int>& result){
                result.clear();                
                if(!current_subset) {
                    for(int at=0;at<Natoms;++at)
                        if( comparison(op1(at),op2(at)) ) result.push_back(at);
                } else {
                    for(int at: *current_subset)
                        if( comparison(op1(at),op2(at)) ) result.push_back(at);
                }                
            };
        }
    }

    else if(node->name == "STR_KEYWORD_EXPR" || node->name == "INT_KEYWORD_EXPR")
    {
        auto keyword = node->nodes[0]->token;

        // Vector of functions used to evaluate atom
        vector<std::function<bool(int)>> evaluators(node->nodes.size()-1);

        // Set evaluators based on the type of value
        for(int i=1;i<node->nodes.size();++i){
            string val = node->nodes[i]->token;
            std::function<bool(int)> func;

            if(node->nodes[i]->name == "STR"){
                // We need direct comparison
                if(keyword=="name")    func = [val,this](int at){ return sys->atoms[at].name == val; };
                else if(keyword=="resname") func = [val,this](int at){ return sys->atoms[at].resname == val; };
                else if(keyword=="tag")     func = [val,this](int at){ return sys->atoms[at].tag == val; };
                else if(keyword=="chain")   func = [val,this](int at){ return sys->atoms[at].chain == val[0]; };
                else if(keyword=="type")    func = [val,this](int at){ return sys->atoms[at].type_name == val; };
            } else if(node->nodes[i]->name == "REGEX") {
                std::regex reg(val);

                if(keyword=="name")    func = [reg,this](int at){ return std::regex_match(sys->atoms[at].name.c_str(),reg); };
                else if(keyword=="resname") func = [reg,this](int at){ return std::regex_match(sys->atoms[at].resname.c_str(),reg); };
                else if(keyword=="tag")     func = [reg,this](int at){ return std::regex_match(sys->atoms[at].tag.c_str(),reg); };
                else if(keyword=="chain")   func = [reg,this](int at){
                    string s(" ");
                    s[0] = sys->atoms[at].chain;
                    return std::regex_match(s.c_str(),reg);
                };
                else if(keyword=="type")    func = [reg,this](int at){ return std::regex_match(sys->atoms[at].type_name.c_str(),reg); };
            } else if(node->nodes[i]->name == "INTEGER") {
                int intval = stol(val);

                if(keyword=="index")    func = [intval](int at){ return at == intval; };
                else if(keyword=="resindex")  func = [intval,this](int at){ return sys->atoms[at].resindex == intval; };
                else if(keyword=="resid")     func = [intval,this](int at){ return sys->atoms[at].resid == intval; };
            } else if(node->nodes[i]->name == "RANGE") {
                int b = stol(node->nodes[i]->nodes[0]->token);
                int e = stol(node->nodes[i]->nodes[1]->token);

                if(keyword=="index")
                    func = [b,e](int at){
                        for(int k=b;k<=e;++k)
                            if(at == k) return true;
                        return false;
                    };
                else if(keyword=="resindex")
                    func = [b,e,this](int at){
                        for(int k=b;k<=e;++k)
                            if(sys->atoms[at].resindex == k) return true;
                        return false;
                    };
                else if(keyword=="resid")
                    func = [b,e,this](int at){
                        for(int k=b;k<=e;++k)
                            if(sys->atoms[at].resid == k) return true;
                        return false;
                    };
            }

            evaluators[i-1] = func;
        }

        returned_function = [evaluators,this](vector<int>& result){
            result.clear();

            if(current_subset){
                for(int at: *current_subset){ // over atoms
                    for(const auto& ev: evaluators) // over values of keyword
                        if( ev(at) ){
                            result.push_back(at);
                            break;
                        }
                }
            } else {
                for(int at=0;at<Natoms;++at){ // over atoms
                    for(const auto& ev: evaluators) // over values of keyword
                        if( ev(at) ){
                            result.push_back(at);
                            break;
                        }
                }
            }            
        };
    }

    else if(node->name == "LOGICAL_EXPR")
    {

        if(node->nodes[1]->token == "or") {
            auto func1 = eval_node(node->nodes[0]);
            auto func2 = eval_node(node->nodes[2]);

            return [func1,func2,this](vector<int>& result){
                result.clear();
                vector<int> res1,res2;
                func1(res1);
                func2(res2);
                std::set_union(res1.begin(),res1.end(),res2.begin(),res2.end(),back_inserter(result));                
            };

        } else if(node->nodes[1]->token == "and") {
            // Optimize to put pure node first
            bool pure1 = !node->nodes[0]->is_coord_dependent;
            bool pure2 = !node->nodes[2]->is_coord_dependent;

            result_func_t func1, func2;
            if(pure2 && !pure1){
                func1 = eval_node(node->nodes[2]);
                func2 = eval_node(node->nodes[0]);
            } else {
                func1 = eval_node(node->nodes[0]);
                func2 = eval_node(node->nodes[2]);
            }

            returned_function = [func1,func2,this](vector<int>& result){
                result.clear();

                vector<int> sub;
                func1(sub);

                vector<int> res2;
                current_subset = &sub; // Set subset
                func2(res2); // Is using filled current subset

                std::set_intersection(sub.begin(),sub.end(),res2.begin(),res2.end(),back_inserter(result));

                // Reset subset
                if(starting_subset){
                    current_subset = starting_subset;
                } else {
                    current_subset = nullptr;
                }                
            };
        }
    }

    else if(node->name == "LOGICAL_OPERAND")
    {
        if(node->nodes[0]->name == "NOT"){
            auto func = eval_node(node->nodes[1]);
            returned_function =  [this,func](vector<int>& result){
                result.clear();
                vector<int> res;
                func(res);
                if(!current_subset){
                    auto r = boost::counting_range(0,Natoms);
                    std::set_difference(r.begin(),r.end(), res.begin(),res.end(), back_inserter(result));
                } else {
                    // For subset
                    std::set_difference(current_subset->begin(),current_subset->end(), res.begin(),res.end(), back_inserter(result));
                }                
            };
        } else if(node->nodes[0]->name == "BYRES"){
            auto func = eval_node(node->nodes[1]);

            returned_function =  [this,func](vector<int>& result){
                result.clear();
                vector<int> res;
                func(res);

                // First make a set of resids we need to search
                std::unordered_set<int> resind;
                for(auto at: res) resind.insert(sys->atoms[at].resindex);

                // Now cycle over all atoms in the starting subset if present (not current subset!!!)
                if(starting_subset){
                    for(int at: *starting_subset) // over starting subset
                        if(resind.count(sys->atoms[at].resindex)) result.push_back(at);
                } else {
                    for(int at=0;at<Natoms;++at) // over all atoms
                        if(resind.count(sys->atoms[at].resindex)) result.push_back(at);
                }
            };
        }
    }

    else if(node->name == "ALL")
    {
        vector<int> res(Natoms);
        for(int at=0;at<Natoms;++at) res[at] = at;

        returned_function =  [res](vector<int>& result){
            result = res;
        };
    }

    else if(node->name == "WITHIN")
    {
        float cutoff = stof(node->nodes[0]->token);

        bool periodic = false;
        bool include_self = true;
        result_func_t func;

        if(node->nodes.size() == 4){ // Both pbc and self
            if(node->nodes[1]->name == "PBC" && (node->nodes[1]->token == "pbc" || node->nodes[1]->token == "periodic"))
                periodic = true;

            if(node->nodes[2]->name == "PBC" && (node->nodes[2]->token == "pbc" || node->nodes[2]->token == "periodic"))
                periodic = true;

            if(node->nodes[1]->name == "SELF" && node->nodes[1]->token == "noself")
                include_self = false;

            if(node->nodes[2]->name == "SELF" && node->nodes[2]->token == "noself")
                include_self = false;

            func = eval_node(node->nodes[3]);

        } else if(node->nodes.size() == 3){ // Either pbc or self
            if(node->nodes[1]->name == "PBC" && (node->nodes[1]->token == "pbc" || node->nodes[1]->token == "periodic"))
                periodic = true;

            if(node->nodes[1]->name == "SELF" && node->nodes[1]->token == "noself")
                include_self = false;

            func = eval_node(node->nodes[2]);

        } else { // Neither pbc nor self
            func = eval_node(node->nodes[1]);
        }

        returned_function =  [this,func,cutoff,include_self,periodic](vector<int>& result){

            Selection dum1(*sys), dum2(*sys);
            // Result is returned directly into the index array of selection dum2
            // thus no additional copying
            func(dum2._index);

            // Prepare selection dum1
            if(!current_subset){
                // We are NOT limited by subspace
                dum1._index.resize(Natoms);
                for(int i=0;i<Natoms;++i) dum1._index[i] = i;
            } else {
                // We are limited by subspace
                dum1._index = *current_subset;
            }

            // Set frame for both selections
            dum1.set_frame(frame);
            dum2.set_frame(frame);

            search_within(cutoff,dum1,dum2,result,include_self,periodic);
        };

    }

    else
    {
        returned_function = [](vector<int>& result){
            result.clear();
        };
    }

    if(!node->is_coord_dependent){
        vector<int> res;
        returned_function(res); // Apply function of itself

        //fmt::print("Precomputing node {}\n",node->name);

        // Return precomputed instead
        returned_function = [res](vector<int>& result){
            result = res;
        };
    }

    return returned_function;

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
std::function<float(int)> Selection_parser::get_numeric(const std::shared_ptr<MyAst> &node){

    std::function<float(int)> res;

    // terminals
    if(node->name == "INTEGER"){
        float val = stol(node->token);
        res = [val](int at){ return val; };
    } else if(node->name == "FLOAT"){
        float val = stof(node->token);
        res = [val](int at){ return val; };
    } else if(node->name == "X"){
        res =[this](int at){ return sys->traj[frame].coord[at](0); };
    } else if(node->name == "Y"){
        res = [this](int at){ return sys->traj[frame].coord[at](1); };
    } else if(node->name == "Z"){
        res = [this](int at){ return sys->traj[frame].coord[at](2); };
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
        auto func = get_numeric(node->nodes[0]);
        res = [func](int at){ return -func(at); };

    } else if(node->name == "NUM_EXPR" || node->name == "NUM_TERM" || node->name == "NUM_POWER") {

        auto func1 = get_numeric(node->nodes[0]);
        auto op = node->nodes[1]->token;
        auto func2 = get_numeric(node->nodes[2]);

        // Evaluation function
        std::function<float(float,float)> func;

        if     (op=="+") {
            func = [](float a,float b){ return a+b; };
        } else if(op=="-") {
            func = [](float a,float b){ return a-b; };
        } else if(op=="*") {
            func = [](float a,float b){ return a*b; };
            res = [func1,func2](int at){ return func1(at)*func2(at); };
        } else if(op=="/") {
            func = [](float a,float b){
                if(b==0.0) throw Pteros_error("Division by zero in selection!");
                return a/b;
            };
        } else if(op=="^" || op=="**") {
            func = [](float a,float b){ return std::pow(a,b); };
        }

        if(node->nodes[0]->is_coord_dependent || node->nodes[2]->is_coord_dependent){
            // For non-pure return evaluation function
            res = [func,func1,func2](int at){ return func(func1(at),func2(at)); };
        } else {
            // For pure return result precomputed for atom 0
            float val = func(func1(0),func2(0));
            res = [val](int at){ return val; };
        }

    } else if(node->name == "POINT") {
        Eigen::Vector3f p;

        bool pbc = false;

        if(node->nodes.size()==3){
            p(0) = get_numeric(node->nodes[0])(0);
            p(1) = get_numeric(node->nodes[1])(0);
            p(2) = get_numeric(node->nodes[2])(0);
        } else {
            // With pbc
            if(node->nodes[0]->token == "pbc" || node->nodes[0]->token == "periodic") pbc = true;
            p(0) = get_numeric(node->nodes[1])(0);
            p(1) = get_numeric(node->nodes[2])(0);
            p(2) = get_numeric(node->nodes[3])(0);
        }

        // Return distance
        if(pbc){
            res = [this,p](int at){
                return sys->box(frame).distance(p, sys->traj[frame].coord[at]);
            };
        } else {
            res = [this,&p](int at){
                return (p - sys->traj[frame].coord[at]).norm();
            };
        }

    } else if(node->name == "VECTOR" || node->name == "PLANE") {
        Eigen::Vector3f p,dir;

        bool pbc = false;

        if(node->nodes.size()==6){
            p(0) = get_numeric(node->nodes[0])(0);
            p(1) = get_numeric(node->nodes[1])(0);
            p(2) = get_numeric(node->nodes[2])(0);
            dir(0) = get_numeric(node->nodes[3])(0);
            dir(1) = get_numeric(node->nodes[4])(0);
            dir(2) = get_numeric(node->nodes[5])(0);
        } else {
            // With pbc
            if(node->nodes[0]->token == "pbc" || node->nodes[0]->token == "periodic") pbc = true;
            p(0) = get_numeric(node->nodes[1])(0);
            p(1) = get_numeric(node->nodes[2])(0);
            p(2) = get_numeric(node->nodes[3])(0);
            dir(0) = get_numeric(node->nodes[4])(0);
            dir(1) = get_numeric(node->nodes[5])(0);
            dir(2) = get_numeric(node->nodes[6])(0);
        }

        if(node->name == "VECTOR"){
            // For vector
            res = [this,p,dir,pbc](int at){
                Eigen::Vector3f atom = sys->traj[frame].coord[at];
                // Get vector from p to current atom
                Eigen::Vector3f v = atom - p;
                // Project v onto dir
                v = (v.dot(dir)/dir.squaredNorm())*dir;
                // Get the end point of projection
                v += p;
                // Return distance between atom and v
                if(pbc){
                    return sys->box(frame).distance(atom, v);
                } else {
                    return (atom-v).norm();
                }
            };
        } else {
            // For plane
            res = [this,p,dir,pbc](int at){
                Eigen::Vector3f atom = sys->traj[frame].coord[at];
                // Get vector from p to current atom
                Eigen::Vector3f v = atom - p;
                // Project v onto dir
                v = (v.dot(dir)/dir.squaredNorm())*dir;
                // Get closest point on a plane to atom
                v = atom-v;
                // Return distance between atom and v
                if(pbc){
                    return sys->box(frame).distance(atom, v);
                } else {
                    return (atom-v).norm();
                }
            };
        }

    }  else {
        throw Pteros_error("Wrong numeric node!");
    }

    return res;    
}


