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


#include "selection_parser.h"
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <Eigen/Core>
#include <unordered_set>
#include <regex>
#include <list>

using namespace std;
using namespace pteros;
using namespace Eigen;


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
    LOGICAL_EXPR        <-  LOGICAL_SEQ(LOGICAL_OPERAND,LOGICAL_OPERATOR)
    LOGICAL_OPERATOR    <-  << 'or' / 'and' >_>
    LOGICAL_OPERAND     <-  KEYWORD_EXPR / NUM_COMP / '(' LOGICAL_EXPR ')' / NOT / BY / WITHIN_EXPR / ALL
    NOT                 <-  'not ' LOGICAL_EXPR     {no_ast_opt}
    ALL                 <-  'all'
    BY                  <-  'by ' BY_PROPERTY LOGICAL_EXPR
    BY_PROPERTY         <-  <'residue' / 'chain' / 'mol'>

    NUM_COMP            <-  NUM_EXPR COMP_OPERATOR NUM_EXPR (COMP_OPERATOR NUM_EXPR)?
    COMP_OPERATOR       <-  < '==' / '<=' / '>=' / '<>' / '!=' / '<' / '>' / '=' >

    NUM_EXPR            <-  NUM_EXPR_SEQ(NUM_OPERAND,NUM_OPERATOR)
    NUM_OPERATOR        <-  < [-+/*^] >
    NUM_OPERAND         <-  FLOAT / NUM_ATOM_PROPERTY / DIST_EXPR / '(' NUM_EXPR ')' / UNARY_MINUS
    UNARY_MINUS         <-  '-' NUM_EXPR     {no_ast_opt}

    FLOAT               <-  < INTEGER ('.' [0-9]+ )? ( ('e' / 'E' ) INTEGER )? >
    INTEGER             <-  < ('-' / '+')? [0-9]+ >

    NUM_ATOM_PROPERTY   <-  X / Y / Z / BETA / OCC / RESID / RESINDEX /
                            INDEX / MASS / CHARGE

    X                   <-  'x ' 'of ' VEC3 / 'x'      {no_ast_opt}
    Y                   <-  'y ' 'of ' VEC3 / 'y'      {no_ast_opt}
    Z                   <-  'z ' 'of ' VEC3 / 'z'      {no_ast_opt}
    BETA                <-  'beta'
    OCC                 <-  'occ'
    RESID               <-  'RESID'
    RESINDEX            <-  'resindex'
    INDEX               <-  'index'
    MASS                <-  'mass'
    CHARGE              <-  'charge'

    KEYWORD_EXPR        <-  STR_KEYWORD_EXPR / INT_KEYWORD_EXPR

    STR_KEYWORD_EXPR    <-  STR_KEYWORD (STR / REGEX)+
    STR_KEYWORD         <-  << 'name' / 'resname' / 'tag' / 'chain' / 'type' >_>
    STR                 <-  !('or'/'and') < [a-zA-Z0-9_]+ >

    INT_KEYWORD_EXPR    <-  INT_KEYWORD (RANGE / INTEGER)+
    INT_KEYWORD         <-  << 'resindex' / 'index' / 'resid' >_>
    RANGE               <-  INTEGER ('-'/'to'/':') INTEGER

    WITHIN_EXPR         <-  'within ' NUM_EXPR (PBC SELF / SELF PBC / PBC / SELF)? 'of ' (LOGICAL_OPERAND / VEC3_EXPR)

    SELF                <-  << 'self' / 'noself' >_>
    PBC                 <-  'pbc ' PBC_DIMS / 'pbc ' / 'nopbc '     {no_ast_opt}
    PBC_DIMS            <-  << [yYnN01]{3} >_>

    DIST_EXPR           <-  'dist ' PBC? 'from ' (VEC3_EXPR / VECTOR / PLANE)       {no_ast_opt}
    VECTOR              <-  'vector ' VEC3_EXPR VEC3_EXPR /
                            'vector ' 'point ' VEC3_EXPR 'dir ' VEC3_EXPR
    PLANE               <-  'plane ' 'point ' VEC3_EXPR 'normal ' VEC3_EXPR /
                            'plane ' VEC3_EXPR VEC3_EXPR VEC3_EXPR

    # Allows vector as 'a b c' or '(a b c)'
    VEC3_EXPR           <-  VEC3 / '(' VEC3 ')'
    VEC3                <-  NUM_EXPR NUM_EXPR NUM_EXPR / CENTER     {no_ast_opt}

    CENTER              <-  'center ' (WEIGHT PBC / PBC WEIGHT / PBC / WEIGHT)? 'of ' LOGICAL_OPERAND {no_ast_opt}
    WEIGHT              <-  'weight ' NUM_EXPR {no_ast_opt}

    %whitespace         <-  [ \t\r\n]*
    _                   <-  [ \t]

    LOGICAL_SEQ(A, O) <-  A (O A)* {
      precedence
        L and or
    }

    NUM_EXPR_SEQ(A, O) <-  A (O A)* {
      precedence
        L + -
        L * /
        R ^
    }

    REGEX               <-  '"' <(!'"' .)*> '"' / "'" <(!"'" .)*> "'"
)");


//===============================================


bool is_node_coordinate_dependent(const std::shared_ptr<MyAst>& node){
    using namespace peg::udl;

    switch(node->tag){
        case "X"_:
        case "Y"_:
        case "Z"_:
        case "WITHIN_EXPR"_:
        case "DIST_EXPR"_:
        case "CENTER"_:
            return true;
        default:
            return false;
    }
}


SelectionParser::SelectionParser(std::vector<int> *subset):
    has_coord(false),
    starting_subset(subset),
    sys(nullptr)
{

}

SelectionParser::~SelectionParser(){}


void set_coord_dependence(const std::shared_ptr<MyAst>& node){
    node->is_coord_dependent = is_node_coordinate_dependent(node);
    if(node->nodes.size()){
        // check children recursively
        for(auto& child: node->nodes){
            set_coord_dependence(child);
            if(child->is_coord_dependent) node->is_coord_dependent = true;
        }
    }
}


void SelectionParser::optimize_numeric(std::shared_ptr<MyAst>& node){
    using namespace peg::udl;

    // optimize arithmetics    
    switch(node->tag) {
    case "NUM_EXPR_SEQ"_:
        // If not-coord dependent just optimize
        if(!node->is_coord_dependent){
            // Replace with float node
            node = std::make_shared<MyAst>("",0,0,"FLOAT", fmt::format("{}", get_numeric(node)(0)));
            node->is_coord_dependent = false; // Keep correct flag just in case
        }
        break;
    }

    // Recurse into children
    for(auto& child: node->nodes) optimize_numeric(child);
}


void SelectionParser::precompute(std::shared_ptr<MyAst>& node){
    using namespace peg::udl;

    switch(node->tag){
    case "NUM_COMP"_:
    case "STR_KEYWORD_EXPR"_:
    case "INT_KEYWORD_EXPR"_:
    case "LOGICAL_SEQ"_:
    case "LOGICAL_OPERAND"_:
    case "ALL"_:
    case "WITHIN_EXPR"_:
    case "BY"_:
    case "NOT"_:
        if(!node->is_coord_dependent){
            //cout << "precomputing " << node->name << endl;
            auto new_node = std::make_shared<MyAst>("",0,0,"PRE","");
            eval_node(node, new_node->precomputed);
            node = new_node;
        } else if(!node->nodes.empty()) {
            // Try children, some of them could be coord-independent
            for(auto& child: node->nodes) precompute(child);
        }
    }
}

void SelectionParser::create_ast(string& sel_str, System* system){
    if (_parser.parse(sel_str.c_str(), tree)) {
        tree = _parser.optimize_ast(tree);
        set_coord_dependence(tree);
    } else {
        throw PterosError(_parser.error_message);
    }

    if(tree->is_coord_dependent) has_coord = true; // Global coord dependence

    sys = system;
    Natoms = sys->num_atoms();

    if(starting_subset)
        current_subset = starting_subset;
    else
        current_subset = nullptr;

    // Optimize numeric values in the tree
    optimize_numeric(tree);

    // proceed with optimizing pure nodes to precomputed if needed
    if(has_coord) precompute(tree);
}

void SelectionParser::apply_ast(size_t fr, vector<int>& result){
    frame = fr;
    eval_node(tree,result);
}


Eigen::Vector3i process_pbc(const std::shared_ptr<MyAst> &node){
    using namespace peg::udl;

    if(node->tag != "PBC"_) throw PterosError("Invalid PBC node!");

    Vector3i ret;
    if(node->choice == 2){ // nopbc
        ret << 0,0,0;
    } else if(node->choice == 1) { // pbc
        ret << 1,1,1;
    } else { // pbc dims
        for(int dim=0;dim<3;++dim){
            char c = node->nodes[0]->token[dim];
            ret(dim) = (c=='y' || c=='Y' || c=='1') ? 1 : 0;
        }
    }
    return ret;
}

void SelectionParser::eval_node(const std::shared_ptr<MyAst> &node, std::vector<int>& result){
    using namespace peg::udl;

    result.clear();

    //==========================
    // Here evaluation starts
    //==========================
    switch(node->tag){
    //---------------------------------------------------------------------------
    case "PRE"_:
        // Precomputed nodes are created during optimization stage
        result = node->precomputed;
        break;

    //---------------------------------------------------------------------------
    case "NUM_COMP"_:
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

        break;
    }

    //---------------------------------------------------------------------------
    case "STR_KEYWORD_EXPR"_:
    {
        std::function<bool(int,const string&)> comp_func_str;
        std::function<bool(int,const std::regex&)> comp_func_regex;

        const auto& keyword = node->nodes[0]->token;

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
                regex_values.emplace_back(string(node->nodes[i]->token));
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

        break;
    }

    //---------------------------------------------------------------------------
    case "INT_KEYWORD_EXPR"_:
    {
        const auto& keyword = node->nodes[0]->token;
        int Nchildren = node->nodes.size(); // Get number of children

        // If starting subset is present than this is a subselection and
        // we have to interpret indexes as local indexes!
        if(keyword == "index") {
            // Cycle over children
            for(int i=1;i<Nchildren;++i){
                if(node->nodes[i]->name == "INTEGER") {
                    int k = node->nodes[i]->token_to_number<int>();
                    // Shift to local index for subselection if needed
                    if(starting_subset) k+=(*starting_subset)[0];
                    // We have to check the range here
                    if(k>=0 && k<Natoms)
                        result.push_back(k);
                } else {
                    // this is a range, not an integer
                    int i1 = node->nodes[i]->nodes[0]->token_to_number<int>();
                    int i2 = node->nodes[i]->nodes[1]->token_to_number<int>();
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
                    int_list.push_back(node->nodes[i]->token_to_number<int>());
                } else {
                    range_list.push_back(node->nodes[i]->nodes[0]->token_to_number<int>());
                    range_list.push_back(node->nodes[i]->nodes[1]->token_to_number<int>());
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

        break;
    }

    //---------------------------------------------------------------------------
    case "LOGICAL_SEQ"_:
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

        break;
    }

    //---------------------------------------------------------------------------
    case "NOT"_:
    {
        vector<int> res;
        eval_node(node->nodes[0],res);

        if(!current_subset){
            vector<int> v(Natoms);
            for(int i=0;i<Natoms;++i) v[i]=i;
            std::set_difference(v.begin(),v.end(), res.begin(),res.end(), back_inserter(result));
        } else {
            // For subset
            std::set_difference(current_subset->begin(),current_subset->end(), res.begin(),res.end(), back_inserter(result));
        }
        break;
    }

    //---------------------------------------------------------------------------
    case "BY"_:
    {
        // Evaluate inner
        vector<int> res;
        eval_node(node->nodes[1],res);

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
            if(!sys->force_field.ready) throw PterosError("Can't select by molecule: no topology!");

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

        break;
    }

    //---------------------------------------------------------------------------
    case "ALL"_:
        result.resize(Natoms);
        for(int at=0;at<Natoms;++at) result[at] = at;

        break;

    //---------------------------------------------------------------------------
    case "WITHIN_EXPR"_:
    {        
        Vector3i pbc = noPBC;
        bool include_self = true;

        // Child 0 is always a cutoff
        // Numeric expression should not be coord dependent!
        if(node->nodes[0]->is_coord_dependent) throw PterosError("Within cutoff can't depend on atomic coordinates!");
        float cutoff = get_numeric(node->nodes[0])(0);

        // Child 1 could be either PBC or SELF
        if(node->nodes[1]->tag == "PBC"_){
            pbc = process_pbc(node->nodes[1]);
            // Check next for self
            if(node->nodes[2]->tag == "SELF"_){
                if(node->nodes[2]->token == "self") include_self = true;
            }
        } else if(node->nodes[1]->tag == "SELF"_) {
            if(node->nodes[1]->token == "self") include_self = true;
            // Check next for pbc
            if(node->nodes[2]->tag == "PBC"_){
                pbc = process_pbc(node->nodes[2]);
            }
        }

        // Prepare first (outer) selection
        Selection dum1(*sys);
        if(!current_subset){
            // We are NOT limited by subspace
            dum1._index.resize(Natoms);
            for(int i=0;i<Natoms;++i) dum1._index[i] = i;
        } else {
            // We are limited by subspace
            dum1._index = *current_subset;
        }
        dum1.set_frame(frame);

        if(node->nodes.back()->tag == "VEC3"_){
            // Distance from point (with abs indexes!)
            DistanceSearchWithin searcher(cutoff,dum1,true,pbc);
            Eigen::Vector3f point = get_vector(node->nodes.back());
            searcher.search_within(point,result);
        } else {
            // Distance between selections
            // Prepare second (inner) selection
            Selection dum2(*sys);
            // Result is returned directly into the index array of selection dum2
            // thus no additional copying
            eval_node(node->nodes.back(),dum2._index);
            dum2.set_frame(frame);
            search_within(cutoff,dum1,dum2,result,include_self,pbc);
        }

        break;
    }

    //---------------------------------------------------------------------------
    default:
        throw PterosError("Unknown node {}!",node->name);   
    } // case
}

//returns a 3-vector
Eigen::Vector3f SelectionParser::get_vector(const std::shared_ptr<MyAst> &node)
{
    using namespace peg::udl;

    switch(node->tag){
    //---------------------------------------------------------------------------
    case "VEC3"_: {
        if(node->nodes.size() == 3){
            // Get 3 floats
            Eigen::Vector3f v;
            v(0) = node->nodes[0]->token_to_number<float>();
            v(1) = node->nodes[1]->token_to_number<float>();
            v(2) = node->nodes[2]->token_to_number<float>();
            return v;
        } else if(node->nodes.size() == 1) {
            // Recurse inside
            return get_vector(node->nodes[0]);
        }
    }

    //---------------------------------------------------------------------------
    case "CENTER"_: {
        auto pbc = noPBC;
        std::function<float(int)> weight_func;

        // Child 0 could be either PBC or WEIGHT
        if(node->nodes[0]->tag == "PBC"_){
            pbc = process_pbc(node->nodes[0]);
            // Check next for weight
            if(node->nodes[1]->tag == "WEIGHT"_){
                weight_func = get_numeric(node->nodes[1]);
            }
        } else if(node->nodes[0]->tag == "WEIGHT"_) {
            weight_func = get_numeric(node->nodes[0]);
            // Check next for pbc
            if(node->nodes[1]->tag == "PBC"_){
                pbc = process_pbc(node->nodes[1]);
            }
        }

        // Create selection to get the center of
        Selection sel(*sys);
        // We have to ignore current subset here!
        // Reset subset
        auto old_subset = current_subset;
        current_subset = nullptr;
        // Evaluate last node directly to selection
        eval_node(node->nodes.back(), sel._index);
        // Set current frame for selection
        sel.set_frame(frame);

        Eigen::Vector3f res;
        if(weight_func){
            // Create vector of weights if nedded
            vector<float> w(sel.size());
            for(int i=0;i<sel.size();++i) w[i] = weight_func(sel.index(i));
            res = sel.center(w,pbc);
        } else {
            // N0 weights
            res = sel.center(false,pbc);
        }

        // Reset subset
        current_subset = old_subset;
        return res;
    }

    //---------------------------------------------------------------------------
    default:
        throw PterosError("Unknown node {}!",node->name);
    } //case

}

// Returns callable, which returns value for numeric node for atom at
std::function<float(int)> SelectionParser::get_numeric(const std::shared_ptr<MyAst> &node){
    using namespace peg::udl;

    switch(node->tag){
    //---------------------------------------------------------------------------
    case "INTEGER"_:
    {
        float val = stol(string(node->token));
        return [val](int at){ return val; };
    }
    //---------------------------------------------------------------------------
    case "FLOAT"_:
    {
        float val = stof(string(node->token));
        return [val](int at){ return val; };
    }
    //---------------------------------------------------------------------------
    case "X"_:
        if(node->nodes.empty())
            return [this](int at){ return sys->traj[frame].coord[at](0); };
        else {
            float x = get_vector(node->nodes[0])[0];
            return [x](int at){ return x; };
        }
    //---------------------------------------------------------------------------
    case "Y"_:
        if(node->nodes.empty())
            return [this](int at){ return sys->traj[frame].coord[at](1); };
        else {
            float y = get_vector(node->nodes[0])[1];
            return [y](int at){ return y; };
        }
    //---------------------------------------------------------------------------
    case "Z"_:
        if(node->nodes.empty())
            return [this](int at){ return sys->traj[frame].coord[at](2); };
        else {
            float z = get_vector(node->nodes[0])[2];
            return [z](int at){ return z; };
        }
    //---------------------------------------------------------------------------
    case "BETA"_:
        return [this](int at){ return sys->atoms[at].beta; };
    //---------------------------------------------------------------------------
    case "OCC"_:
        return [this](int at){ return sys->atoms[at].occupancy; };
    //---------------------------------------------------------------------------
    case "INDEX"_:
        return [](int at){ return at; };
    //---------------------------------------------------------------------------
    case "RESINDEX"_:
        return [this](int at){ return sys->atoms[at].resindex; };
    //---------------------------------------------------------------------------
    case "RESID"_:
        return [this](int at){ return sys->atoms[at].resid; };
    //---------------------------------------------------------------------------
    case "MASS"_:
        return [this](int at){ return sys->atoms[at].mass; };
    //---------------------------------------------------------------------------
    case "CHARGE"_:
        return [this](int at){ return sys->atoms[at].charge; };
    //---------------------------------------------------------------------------
    // Compounds
    case "UNARY_MINUS"_:
    {
        auto func = get_numeric(node->nodes[0]);
        return [func](int at){ return -func(at); };
    }
    //---------------------------------------------------------------------------
    case "NUM_EXPR_SEQ"_:
    {        
        std::function<float(int)> operand1,operand2;
        std::function<float(float,float)> op;

        operand1 = get_numeric(node->nodes[0]);
        operand2 = get_numeric(node->nodes[2]);

        if(node->nodes[1]->token == "+")
            op= [](float a,float b){ return a+b; };
        else if(node->nodes[1]->token == "-")
            op= [](float a,float b){ return a-b; };
        else if(node->nodes[1]->token == "*")
            op= [](float a,float b){ return a*b; };
        else if(node->nodes[1]->token == "/")
            op= [](float a,float b){
                if(b==0.0) throw PterosError("Division by zero in selection!");
                return a/b;
            };            
        else if(node->nodes[1]->token == "^")
            op = [](float a,float b){ return std::pow(a,b); };

        return [operand1,operand2,op](int at){
            return op(operand1(at),operand2(at));
        };
    }
    //---------------------------------------------------------------------------
    case "DIST_EXPR"_:
    {
        Array3i pbc = noPBC;
        int offset = 0;

        // Process PBC if present
        if(node->nodes[0]->tag == "PBC"_){
            pbc = process_pbc(node->nodes[0]);
            ++offset;
        }

        // Inner node
        const auto& inner = node->nodes[offset];

        switch(inner->tag){

        case "VEC3"_: { // from point
            Eigen::Vector3f p = get_vector(inner);

            // Return distance
            if((pbc!=0).any()){
                return [this,p](int at){
                    return sys->box(frame).distance(p, sys->traj[frame].coord[at]);
                };
            } else {
                return [this,p](int at){
                    return (p - sys->traj[frame].coord[at]).norm();
                };
            }

            break;
        }

        case "VECTOR"_: { // from vector
            Eigen::Vector3f p,dir;
            if(inner->choice == 0){ // Two points
                p = get_vector(inner->nodes[0]);
                dir = (get_vector(inner->nodes[1]) - p).normalized();
            } else { // point and direction
                p = get_vector(inner->nodes[0]);
                dir = get_vector(inner->nodes[1]).normalized();
            }

            if((pbc.array()!=0).any()){ // periodic
                return [this,p,dir,pbc](int at){
                    Eigen::Vector3f atom = sys->traj[frame].coord[at];
                    // Get vector from p to current atom
                    Eigen::Vector3f v = atom - p;
                    // Project v onto dir
                    v = (v.dot(dir)/dir.squaredNorm())*dir;
                    // Get the end point of projection
                    v += p;
                    // Return periodic distance between atom and v
                    return sys->box(frame).distance(atom, v, pbc);
                };
            } else { // nonperiodic
                return [this,p,dir,pbc](int at){
                    Eigen::Vector3f atom = sys->traj[frame].coord[at];
                    // Get vector from p to current atom
                    Eigen::Vector3f v = atom - p;
                    // Project v onto dir
                    v = (v.dot(dir)/dir.squaredNorm())*dir;
                    // Get the end point of projection
                    v += p;
                    // Return distance between atom and v
                    return (atom-v).norm();
                };
            }

            break;
        }

        case "PLANE"_: { // From plane
            Eigen::Vector3f p,dir;
            if(inner->choice == 0){ // point and normal
                p = get_vector(inner->nodes[0]);
                dir = get_vector(inner->nodes[1]).normalized();
            } else { // 3 points
                p = get_vector(inner->nodes[0]);
                Eigen::Vector3f vec1 = get_vector(inner->nodes[1])-p;
                Eigen::Vector3f vec2 = get_vector(inner->nodes[2])-p;
                dir = vec1.cross(vec2).normalized();
            }

            if((pbc.array()!=0).any()){ // periodic
                return [this,p,dir,pbc](int at){
                    Eigen::Vector3f atom = sys->traj[frame].coord[at];
                    // Get vector from p to current atom
                    Eigen::Vector3f v = atom - p;
                    // Project v onto dir
                    v = (v.dot(dir)/dir.squaredNorm())*dir;
                    // Get closest point on a plane to atom
                    v = atom-v;
                    // Return distance between atom and v
                    return sys->box(frame).distance(atom, v, pbc);
                };
            } else { // nonperiodic
                return [this,p,dir](int at){
                    Eigen::Vector3f atom = sys->traj[frame].coord[at];
                    // Get vector from p to current atom
                    Eigen::Vector3f v = atom - p;
                    // Project v onto dir
                    v = (v.dot(dir)/dir.squaredNorm())*dir;
                    // Get closest point on a plane to atom
                    v = atom-v;
                    // Return distance between atom and v
                    return (atom-v).norm();
                };
            }

        }
        } // switch


    } // DIST_EXPR
    //---------------------------------------------------------------------------
    default:
        throw PterosError("Wrong numeric node!");

    } // case
}





