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
#include <boost/algorithm/string.hpp>
#include <boost/range/counting_range.hpp>
#include <unordered_set>
#include <regex>
#include <list>

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

        X                  <-  ('X' / 'x') ('of' COM)?
        Y                  <-  ('Y' / 'y') ('of' COM)?
        Z                  <-  ('Z' / 'z') ('of' COM)?
        BETA               <- < 'beta' >
        OCC                <- < 'occupancy' / 'occ' >
        RESINDEX           <- < 'resindex' >
        INDEX              <- < 'index' >
        RESID              <- < 'resid' >

        COM                <- COM_TYPE PBC? 'of' LOGICAL_OPERAND
        COM_TYPE           <- < 'com' / 'cog' >

        DIST               <- ('dist' / 'distance') (POINT / VECTOR / PLANE)
        POINT              <- 'point' PBC? (FLOAT FLOAT FLOAT / COM)
        VECTOR             <- 'vector' PBC? (FLOAT FLOAT FLOAT / COM) FLOAT FLOAT FLOAT
        PLANE              <- 'plane' PBC? (FLOAT FLOAT FLOAT / COM) FLOAT FLOAT FLOAT

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
            || node->name == "POINT" || node->name == "PLANE" || node->name == "VECTOR" || node->name == "COM"
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
    if(node->nodes.size()){
        // tree
        node->is_coord_dependent = is_node_coordinate_dependent(node);
        for(int i=0;i<node->nodes.size();++i){
            set_coord_dependence(node->nodes[i]);
            if(node->nodes[i]->is_coord_dependent) node->is_coord_dependent = true;
        }
    }

    //cout << node->name << " :: " << node->is_coord_dependent << endl;
}

void Selection_parser::precompute(std::shared_ptr<MyAst>& node){
    if(    node->name!="PRE"
        && node->name!="NUM_COMPARISON"
        && node->name!="STR_KEYWORD_EXPR"
        && node->name!="INT_KEYWORD_EXPR"
        && node->name!="LOGICAL_EXPR"
        && node->name!="LOGICAL_OPERAND"
        && node->name!="ALL"
        && node->name!="WITHIN"
        && node->name!="BYRES") return;

    if(!node->is_coord_dependent){
        node->nodes.clear();
        auto ast = std::make_shared<MyAst>(*node,"PRE");
        eval_node(ast, ast->precomputed);
        node = ast;
    } else {
        for(int i=0;i<node->nodes.size();++i) precompute(node->nodes[i]);
    }
}

void Selection_parser::create_ast(string& sel_str, System* system){
    if (_parser.parse(sel_str.c_str(), tree)) {
        tree = peg::AstOptimizer(true,{"POINT","X","Y","Z"}).optimize(tree);

        //cout << peg::ast_to_s(tree);

        set_coord_dependence(tree);
    } else {
        throw Pteros_error(_parser.error_message);
    }

    if(tree->is_coord_dependent) has_coord = true; // Global coord dependence

    sys = system;
    Natoms = sys->num_atoms();

    if(starting_subset)
        current_subset = starting_subset;
    else
        current_subset = nullptr;

    // proceed with optimizing pure nodes to precomputed if needed
    if(has_coord) precompute(tree);
}

void Selection_parser::apply_ast(size_t fr, vector<int>& result){
    frame = fr;
    eval_node(tree,result);
}


void Selection_parser::eval_node(const std::shared_ptr<MyAst> &node, std::vector<int>& result){

    result.clear();

    // Here starts evaluation

    if(node->name == "PRE"){ // Precomputed nodes are created during optimization stage
        result = node->precomputed;
    }

    else if(node->name == "NUM_COMPARISON")
    {
        auto op  = node->nodes[1]->token;
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
            result.clear();

        } else {
            result.clear();
            if(!current_subset) {
                for(int at=0;at<Natoms;++at)
                    if( comparison(op1(at),op2(at)) ) result.push_back(at);
            } else {
                for(int at: *current_subset)
                    if( comparison(op1(at),op2(at)) ) result.push_back(at);
            }
        }
    }

    else if(node->name == "STR_KEYWORD_EXPR")
    {
        std::function<bool(int,const string&)> comp_func_str;
        std::function<bool(int,const std::regex&)> comp_func_regex;

        const string& keyword = node->nodes[0]->token;

        if(keyword == "name"){
            comp_func_str = [this](int at, const string& str){ return sys->atoms[at].name == str; };
            comp_func_regex = [this](int at, const std::regex& reg){ return std::regex_match(sys->atoms[at].name.c_str(),reg); };
        } else if(keyword == "type"){
            comp_func_str = [this](int at, const string& str){ return sys->atoms[at].type_name == str; };
            comp_func_regex = [this](int at, const std::regex& reg){ return std::regex_match(sys->atoms[at].type_name.c_str(),reg); };
        } else if(keyword == "resname"){
            comp_func_str = [this](int at, const string& str){ return sys->atoms[at].resname == str; };
            comp_func_regex = [this](int at, const std::regex& reg){ return std::regex_match(sys->atoms[at].resname.c_str(),reg); };
        } else if(keyword == "tag"){
            comp_func_str = [this](int at, const string& str){ return sys->atoms[at].tag== str; };
            comp_func_regex = [this](int at, const std::regex& reg){ return std::regex_match(sys->atoms[at].tag.c_str(),reg); };
        } else if(keyword == "chain"){
            comp_func_str = [this](int at, const string& str){ return sys->atoms[at].chain == str[0]; };
            comp_func_regex = [this](int at, const std::regex& reg){
                string s(" ");
                s[0] = sys->atoms[at].chain;
                return std::regex_match(s.c_str(),reg);
            };
        }

        list<string> str_values;
        list<std::regex> regex_values;
        for(int i=1;i<node->nodes.size();++i){
            if(node->nodes[i]->name=="STR")
                str_values.emplace_back(node->nodes[i]->token);
            else
                regex_values.emplace_back(node->nodes[i]->token);
        }

        // Loop body
        auto body = [&](int at){
            // Cycle over regex values
            bool reg_found = false;
            for(const auto& reg: regex_values){
                if(comp_func_regex(at,reg)){
                    result.push_back(at);
                    reg_found = true;
                    break;
                }
            }

            // If at least one regex matched no need to proceed with strings
            if(reg_found) return;

            // Cycle over string values
            for(const auto& str: str_values){
                if(comp_func_str(at,str)){
                    result.push_back(at);
                    break;
                }
            }
        };

        if(!current_subset){
            for(int at=0;at<Natoms;++at) body(at);
        } else {
            for(int at: *current_subset) body(at);
        }

        // Sort unique
        sort(result.begin(),result.end());
        vector<int>::iterator it = unique(result.begin(), result.end());
        result.resize( it - result.begin() );
    }

    //---------------------------------------------------------------------------
    else if(node->name == "INT_KEYWORD_EXPR")
    {
        const string& keyword = node->nodes[0]->token;
        int Nchildren = node->nodes.size(); // Get number of children

        if(keyword == "index") {
            // Cycle over children
            for(int i=0;i<Nchildren;++i){
                if(node->nodes[i]->name == "integer") {
                    int k = stoi(node->nodes[i]->token);
                    // We have to check the range here
                    if(k>=0 && k<Natoms)
                        result.push_back(k);
                } else {
                    // this is a range, not an integer
                    int i1 = stoi(node->nodes[i]->nodes[0]->token);
                    int i2 = stoi(node->nodes[i]->nodes[1]->token);
                    for(int k=i1;k<=i2;++k)
                        // We have to check the range here
                        if(k>=0 && k<Natoms)
                            result.push_back(k);
                }
            }

        } else { // resid or resindex

            // Make lists
            vector<int> int_list;
            vector<int> range_list;
            for(int i=1;i<Nchildren;++i){
                if(node->nodes[i]->name == "INTEGER") {
                    int_list.push_back(stoi(node->nodes[i]->token));
                } else {
                    range_list.push_back(stoi(node->nodes[i]->nodes[0]->token));
                    range_list.push_back(stoi(node->nodes[i]->nodes[1]->token));
                }
            }

            // Comparison function
            std::function<bool(int,int)> comp_func;
            if(keyword == "resid"){
                comp_func = [this](int at, int k){ return sys->atoms[at].resid == k; };
            } else if(keyword == "resindex"){
                comp_func = [this](int at, int k){ return sys->atoms[at].resindex == k; };
            }

            // Loop body
            auto body = [&](int at){
                // Individual numbers
                for(int k: int_list)
                    if(comp_func(at,k)) result.push_back(at);
                // Ranges
                for(int i=0;i<range_list.size();i+=2){ // Itarage by pair
                    for(int k=range_list[i];k<=range_list[i+1];++k){ // Inside range
                        if(comp_func(at,k)) result.push_back(at);
                    }
                }
            };

            // Do loop
            if(!current_subset){
                for(int at=0;at<Natoms;++at) body(at);
            } else {
                for(int at: *current_subset) body(at);
            }
        } // index if

        sort(result.begin(),result.end());
        vector<int>::iterator it = unique(result.begin(), result.end());
        result.resize( it - result.begin() );
    }

    //---------------------------------------------------------------------------

    else if(node->name == "LOGICAL_EXPR")
    {
        if(node->nodes[1]->token == "or") {
            vector<int> res1,res2;
            eval_node(node->nodes[0],res1);
            eval_node(node->nodes[2],res2);
            std::set_union(res1.begin(),res1.end(),res2.begin(),res2.end(),back_inserter(result));

        } else if(node->nodes[1]->token == "and") {
            // Optimize to put pure node first
            bool pure1 = !node->nodes[0]->is_coord_dependent;
            bool pure2 = !node->nodes[2]->is_coord_dependent;

            if(pure2 && !pure1) std::swap(node->nodes[0],node->nodes[2]);

            vector<int> res1,res2;
            eval_node(node->nodes[0],res1);
            current_subset = &res1; // Set subset for second
            eval_node(node->nodes[2],res2); // Is using filled current subset

            std::set_intersection(res1.begin(),res1.end(),res2.begin(),res2.end(),back_inserter(result));

            // Reset subset
            if(starting_subset){
                current_subset = starting_subset;
            } else {
                current_subset = nullptr;
            }
        }
    }

    else if(node->name == "LOGICAL_OPERAND")
    {
        vector<int> res;
        eval_node(node->nodes[1],res);

        if(node->nodes[0]->name == "NOT"){
            if(!current_subset){
                auto r = boost::counting_range(0,Natoms);
                std::set_difference(r.begin(),r.end(), res.begin(),res.end(), back_inserter(result));
            } else {
                // For subset
                std::set_difference(current_subset->begin(),current_subset->end(), res.begin(),res.end(), back_inserter(result));
            }

        } else if(node->nodes[0]->name == "BYRES"){
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
        }
    }

    else if(node->name == "ALL")
    {
        result.resize(Natoms);
        for(int at=0;at<Natoms;++at) result[at] = at;
    }

    else if(node->name == "WITHIN")
    {        
        float cutoff = stof(node->nodes[0]->token);
        bool periodic = false;
        bool include_self = true;        

        int eval_ind;

        if(node->nodes.size() == 4){ // Both pbc and self
            if(node->nodes[1]->name == "PBC" && (node->nodes[1]->token == "pbc" || node->nodes[1]->token == "periodic"))
                periodic = true;

            if(node->nodes[2]->name == "PBC" && (node->nodes[2]->token == "pbc" || node->nodes[2]->token == "periodic"))
                periodic = true;

            if(node->nodes[1]->name == "SELF" && node->nodes[1]->token == "noself")
                include_self = false;

            if(node->nodes[2]->name == "SELF" && node->nodes[2]->token == "noself")
                include_self = false;

            eval_ind = 3;

        } else if(node->nodes.size() == 3){ // Either pbc or self
            if(node->nodes[1]->name == "PBC" && (node->nodes[1]->token == "pbc" || node->nodes[1]->token == "periodic"))
                periodic = true;

            if(node->nodes[1]->name == "SELF" && node->nodes[1]->token == "noself")
                include_self = false;

            eval_ind = 2;

        } else { // Neither pbc nor self
            eval_ind = 1;
        }

        Selection dum1(*sys), dum2(*sys);
        // Result is returned directly into the index array of selection dum2
        // thus no additional copying
        eval_node(node->nodes[eval_ind],dum2._index);

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
    }

    else
    {
        throw Pteros_error("Unknown node {}!",node->name);
    }
}

//returns a 3-vector
Eigen::Vector3f Selection_parser::get_vector(const std::shared_ptr<MyAst> &node){
    if(node->name == "COM")
    {
        bool with_mass = (node->nodes[0]->token == "com") ? true : false;

        auto pbc = noPBC;
        if(node->nodes[1]->name == "PBC" && (node->nodes[1]->token == "pbc" || node->nodes[1]->token == "periodic")) pbc = fullPBC;

        Selection sel(*sys);
        eval_node(node->nodes.back(), sel._index);
        sel.set_frame(frame);

        return sel.center(with_mass,pbc);
    }
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
        if(node->nodes.empty())
            res =[this](int at){ return sys->traj[frame].coord[at](0); };
        else {
            float x = get_vector(node->nodes[0])[0];
            res =[x](int at){ return x; };
        }
    } else if(node->name == "Y"){
        if(node->nodes.empty())
            res =[this](int at){ return sys->traj[frame].coord[at](1); };
        else {
            float y = get_vector(node->nodes[0])[1];
            res =[y](int at){ return y; };
        }
    } else if(node->name == "Z"){
        if(node->nodes.empty())
            res =[this](int at){ return sys->traj[frame].coord[at](2); };
        else {
            float z = get_vector(node->nodes[0])[2];
            res =[z](int at){ return z; };
        }
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

        int N = node->nodes.size();
        int offset = 0;

        bool pbc = false;
        if(node->nodes[0]->token == "pbc" || node->nodes[0]->token == "periodic"){
            pbc = true;
            offset = 1;
        }

        if(N==1 || N==2){ // com and possibly pbc
            p = get_vector(node->nodes[0+offset]);
        } else { // 3 floats and possibly pbc
            p(0) = get_numeric(node->nodes[0+offset])(0);
            p(1) = get_numeric(node->nodes[1+offset])(0);
            p(2) = get_numeric(node->nodes[2+offset])(0);
        }

        // Return distance
        if(pbc){
            res = [this,p](int at){
                return sys->box(frame).distance(p, sys->traj[frame].coord[at]);
            };
        } else {
            res = [this,p](int at){
                return (p - sys->traj[frame].coord[at]).norm();
            };
        }

    } else if(node->name == "VECTOR" || node->name == "PLANE") {
        Eigen::Vector3f p,dir;

        int N = node->nodes.size();
        int offset = 0;

        bool pbc = false;
        if(node->nodes[0]->token == "pbc" || node->nodes[0]->token == "periodic"){
            pbc = true;
            offset = 1;
        }

        if(N==4 || N==5){ // com + 3 floats
            p = get_vector(node->nodes[0+offset]);
            dir(0) = get_numeric(node->nodes[1+offset])(0);
            dir(1) = get_numeric(node->nodes[2+offset])(0);
            dir(2) = get_numeric(node->nodes[3+offset])(0);
        } if(N==6){ // 6 floats
            p(0) = get_numeric(node->nodes[0+offset])(0);
            p(1) = get_numeric(node->nodes[1+offset])(0);
            p(2) = get_numeric(node->nodes[2+offset])(0);
            dir(0) = get_numeric(node->nodes[3+offset])(0);
            dir(1) = get_numeric(node->nodes[4+offset])(0);
            dir(2) = get_numeric(node->nodes[5+offset])(0);
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


