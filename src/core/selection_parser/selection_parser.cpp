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
        LOGICAL_OPERAND    <-  (NOT / BY)? ( '(' LOGICAL_EXPR ')' / ALL / NUM_COMPARISON / KEYWORD_EXPR / WITHIN )
        ALL                <-  < 'all' >
        NOT                <-  < 'not' >
        BY                 <-  'by' < ('residue' / 'chain' / 'mol') > / 'same' < ('residue' / 'chain' / 'mol') > 'as'

        NUM_COMPARISON     <-  NUM_EXPR COMPARISON_OPERATOR NUM_EXPR (COMPARISON_OPERATOR NUM_EXPR)?
        COMPARISON_OPERATOR <- <  '==' / '<=' / '>=' / '<>' / '!=' / '<' / '>' / '=' >
        NUM_EXPR           <- NUM_TERM (PLUS_MINUS NUM_TERM)*
        NUM_TERM           <- NUM_POWER (DIV_MUL NUM_POWER)*
        NUM_POWER          <- NUM_FACTOR (POW NUM_FACTOR)?
        NUM_FACTOR         <- UNARY_MINUS? ( '(' NUM_EXPR ')' / X / Y / Z / BETA / OCC / RESINDEX / INDEX
                                                              / RESID / DIST / MASS / CHARGE ) / FLOAT
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
        MASS               <- < 'mass' >
        CHARGE             <- < 'charge' >

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
        INT_KEYWORD        <- < 'resindex' / 'index' / 'resid' >
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


void Selection_parser::optimize(std::shared_ptr<MyAst>& node){
    // optimize arithmetics
    if(node->name == "NUM_EXPR" || node->name == "NUM_TERM" || node->name == "NUM_POWER") {
        // If not-coord dependent just optimize
        if(!node->is_coord_dependent){
            // Replace with float node
            node = std::make_shared<MyAst>("",0,0,"FLOAT", fmt::format("{}", get_numeric(node)(0)));
            node->is_coord_dependent = false; // Keep correct flag just in case
        }
    }

    // Convert chained logic into a tree
    if(node->name == "LOGICAL_EXPR" && node->nodes.size()>3){
        // second operand and all after it should become new node
        vector<std::shared_ptr<MyAst>> ch;
        copy(node->nodes.begin()+2,node->nodes.end(),back_inserter(ch));
        auto operand2 = std::make_shared<MyAst>("",0,0,"LOGICAL_EXPR",ch);
        // Inherit coord dependence for this new node (important in order not to mislead precomputer)
        operand2->is_coord_dependent = node->is_coord_dependent;
        // keep only 3 nodes
        node->nodes.resize(3);
        // Put new second operand
        node->nodes[2] = operand2;
    }

    // Recurse into children
    for(int i=0;i<node->nodes.size();++i) optimize(node->nodes[i]);
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
        auto ast = std::make_shared<MyAst>("",0,0,"PRE","");
        eval_node(node, ast->precomputed);
        node = ast;
    } else if(!node->nodes.empty()) {
        for(int i=0;i<node->nodes.size();++i) precompute(node->nodes[i]);
    }
}

void Selection_parser::create_ast(string& sel_str, System* system){
    if (_parser.parse(sel_str.c_str(), tree)) {
        tree = peg::AstOptimizer(true,{"POINT","X","Y","Z"}).optimize(tree);

        //cout << sel_str << endl;
        //cout << peg::ast_to_s(tree) << endl;

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

    // Optimize tree
    optimize(tree);

    //cout << peg::ast_to_s(tree) << endl;

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
        vector<string> c; // comparison operators
        vector<std::function<float(int)>> op; // comparison operands
        vector<std::function<bool(float,float)>> comparison; // Function(s) to evaluate

        if(node->nodes.size() == 3){ // simple
            c.resize(1);
            op.resize(2);
            comparison.resize(1);

            c[0]  = node->nodes[1]->token;
            op[0] = get_numeric(node->nodes[0]);
            op[1] = get_numeric(node->nodes[2]);

        } else { // chained
            c.resize(2);
            op.resize(3);
            comparison.resize(2);

            c[0]  = node->nodes[1]->token;
            c[1]  = node->nodes[3]->token;
            op[0] = get_numeric(node->nodes[0]);
            op[1] = get_numeric(node->nodes[2]);
            op[2] = get_numeric(node->nodes[4]);
        }

        for(int i=0;i<c.size();++i){
            if(c[i] == "=" || c[i]== "=="){
                comparison[i] = [](float a, float b){ return a==b; };
            } else if (c[i] == "!=" || c[i]=="<>"){
                comparison[i] = [](float a, float b){ return a!=b; };
            } else if (c[i] == "<"){
                comparison[i] = [](float a, float b){ return a<b; };
            } else if (c[i] == ">"){
                comparison[i] = [](float a, float b){ return a>b; };
            } else if (c[i] == "<="){
                comparison[i] = [](float a, float b){ return a<=b; };
            } else if (c[i] == ">="){
                comparison[i] = [](float a, float b){ return a>=b; };
            }
        }

        if(!current_subset) {
            if(node->nodes.size() == 3){ // simple
                for(int at=0;at<Natoms;++at)
                    if( comparison[0](op[0](at),op[1](at)) ) result.push_back(at);
            } else { // chained
                for(int at=0;at<Natoms;++at)
                    if( comparison[0](op[0](at),op[1](at)) && comparison[1](op[1](at),op[2](at)) ) result.push_back(at);
            }
        } else {
            if(node->nodes.size() == 3){ // simple
                for(int at: *current_subset)
                    if( comparison[0](op[0](at),op[1](at)) ) result.push_back(at);
            } else { // chained
                for(int at: *current_subset)
                    if( comparison[0](op[0](at),op[1](at)) && comparison[1](op[1](at),op[2](at)) ) result.push_back(at);
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

        // If starting subset is present than this is a subselection and
        // we have to interpret indexes as local indexes!
        if(keyword == "index") {
            // Cycle over children
            for(int i=1;i<Nchildren;++i){
                if(node->nodes[i]->name == "INTEGER") {
                    int k = stoi(node->nodes[i]->token);
                    // Shift to local index for subselection if needed
                    if(starting_subset) k+=(*starting_subset)[0];
                    // We have to check the range here
                    if(k>=0 && k<Natoms)
                        result.push_back(k);
                } else {
                    // this is a range, not an integer
                    int i1 = stoi(node->nodes[i]->nodes[0]->token);
                    int i2 = stoi(node->nodes[i]->nodes[1]->token);
                    // Shift to local index for subselection if needed
                    if(starting_subset){
                        i1+=(*starting_subset)[0];
                        i2+=(*starting_subset)[0];
                    }
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

        } else if(node->nodes[0]->name == "BY"){
            if(node->nodes[0]->token == "residue"){
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
                std::sort(result.begin(),result.end());

            } else if(node->nodes[0]->token == "chain") {
                // First make a set of chains we need to search
                std::unordered_set<char> chains;
                for(auto at: res) chains.insert(sys->atoms[at].chain);

                // Now cycle over all atoms in the starting subset if present (not current subset!!!)
                if(starting_subset){
                    for(int at: *starting_subset) // over starting subset
                        if(chains.count(sys->atoms[at].chain)) result.push_back(at);
                } else {
                    for(int at=0;at<Natoms;++at) // over all atoms
                        if(chains.count(sys->atoms[at].chain)) result.push_back(at);
                }
                std::sort(result.begin(),result.end());

            } else if(node->nodes[0]->token == "mol") {
                if(!sys->force_field.ready) throw Pteros_error("Can't select by molecule: no topology!");

                // set of mols we need to include
                std::unordered_set<int> mols;

                // go over res and add mol is we need it
                for(auto at: res){
                    // Cycle over all molecules in the system
                    for(int j=0; j<sys->force_field.molecules.size(); ++j){ // Over molecules
                        if(at>=sys->force_field.molecules[j](0) && at<=sys->force_field.molecules[j](1)){
                            mols.insert(j);
                            break;
                        }
                    }
                }

                // Now add all chosen molecules
                for(int m: mols)
                    for(int at=sys->force_field.molecules[m](0); at<=sys->force_field.molecules[m](1); ++at)
                        result.push_back(at);

                // Sort and restrict to starting subset (!) if needed
                sort(result.begin(),result.end());
                if(starting_subset)
                    std::set_intersection(result.begin(),result.end(),
                                          starting_subset->begin(),starting_subset->end(),
                                          back_inserter(result));

            } // mol
        } //BY
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
        // We have to ignore current subset here!
        // Reset subset
        auto old_subset = current_subset;
        current_subset = nullptr;

        eval_node(node->nodes.back(), sel._index);
        sel.set_frame(frame);

        current_subset = old_subset;

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
    } else if(node->name == "MASS"){
        res = [this](int at){ return sys->atoms[at].mass; };
    } else if(node->name == "CHARGE"){
        res = [this](int at){ return sys->atoms[at].charge; };
    // Compounds
    } else if(node->name == "UNARY_MINUS"){
        auto func = get_numeric(node->nodes[0]);
        res = [func](int at){ return -func(at); };

    } else if(node->name == "NUM_EXPR" || node->name == "NUM_TERM" || node->name == "NUM_POWER") {

        int N = node->nodes.size();

        vector<std::function<float(int)>> operands((N-1)/2+1);
        vector<std::function<float(float,float)>> operators((N-1)/2);

        for(int i=0;i<operands.size();++i){
            operands[i] = get_numeric(node->nodes[i*2]);
        }

        for(int i=0;i<operators.size();++i){
            auto op = node->nodes[i*2+1]->token;
            if     (op=="+") {
                operators[i] = [](float a,float b){ return a+b; };
            } else if(op=="-") {
                operators[i] = [](float a,float b){ return a-b; };
            } else if(op=="*") {
                operators[i] = [](float a,float b){ return a*b; };
            } else if(op=="/") {
                operators[i] = [](float a,float b){
                    if(b==0.0) throw Pteros_error("Division by zero in selection!");
                    return a/b;
                };
            } else if(op=="^" || op=="**") {
                operators[i] = [](float a,float b){ return std::pow(a,b); };
            }
        }

        res = [operands,operators](int at){
            float result = operands[0](at);
            for(int i=0;i<operators.size();++i) result = operators[i](result,operands[i+1](at));
            return result;
        };

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


