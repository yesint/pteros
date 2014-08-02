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

#ifndef SELECTION_GRAMMAR_H
#define SELECTION_GRAMMAR_H

#include "pteros/core/selection_parser.h"
#include <functional>
#include <map>
#include "pteros/core/pteros_error.h"

//===========================================================
using namespace std;
using namespace pteros;


typedef std::function<AstNode_ptr()> result_t;

class Parse_node {
public:
#ifdef _DEBUG_PARSER
    std::string name; // Rule name
#endif
    result_t result; // result of the rule (code returning AstNode_ptr)
    std::vector<std::shared_ptr<Parse_node>> children; // Sub-rules    
};

typedef std::shared_ptr<Parse_node> Parse_node_ptr;

#ifdef _DEBUG_PARSER
#define DEBUG(code) code
#else
#define DEBUG(code)
#endif

#define HEADER(_name) \
    if(_pos_==end) return false; \
    Parse_node_ptr _this_rule_(new Parse_node); \
    _this_rule_->name = #_name; \
    Parse_node_ptr saved_parent = current_parent; \
    current_parent = _this_rule_; \
    static int rule_id = -1; \
    if(rule_id==-1){ /*first call of this rule*/ \
        rule_id = rule_counter; \
        rule_counter++; \
        if(memo.size()<rule_id+1) memo.resize(rule_id+1);\
    } \
    int n = std::distance(beg,_pos_); \
    DEBUG(for(int i=0;i<level;++i) cout <<"  ";) \
    DEBUG(cout << "::: " << _this_rule_->name << " id: " << rule_id << " at: " << n << endl;) \
    /* variables for rule body */\
    std::string::iterator _old_ = _pos_; \
    bool _ok_ = false; \
    /* check memotable */ \
    bool restored = false; \
    if(do_memo && memo[rule_id].count(n)==1){ \
        DEBUG(for(int i=0;i<level;++i) cout <<"  ";) \
        DEBUG(cout << "from memo: " << _this_rule_->name << " at: " << n << endl;) \
        Memo_data& m = memo[rule_id][n]; \
        _this_rule_ = m.tree; \
        _pos_ = m.pos; \
        _ok_ = m.ok; \
        restored = true; \
        num_restored++; \
    } \
    level++; \
    if(max_level<level) max_level = level;\
    if(!restored){ \


// Rule
#define RULE(_name) \
bool _name(bool add_to_tree = true, bool do_memo = true){ \
    bool simplify = false;\
    HEADER(_name) \

#define RULE_REDUCE(_name) \
bool _name(bool add_to_tree = true, bool do_memo = true){ \
    bool simplify = true;\
    HEADER(_name) \

#define END_RULE \
    } /*if not restored*/ \
    if(_ok_){    \
        if(add_to_tree) { \
            DEBUG(for(int i=0;i<level;++i) cout <<"  ";) \
            DEBUG(cout << "Adding node "<<_this_rule_->name << " to parent " << saved_parent->name << endl;) \
            last_success = _pos_;\
            if(simplify && _this_rule_->children.size()==1)\
                saved_parent->children.push_back(_this_rule_->children[0]); \
            else \
                saved_parent->children.push_back(_this_rule_); \
        } \
    } else { \
        _pos_ = _old_; \
        _this_rule_->children.clear(); \
        _this_rule_->result = nullptr; \
    } \
    current_parent = saved_parent; \
    /* Add to memotable */ \
    if(do_memo && !restored) { \
       DEBUG(for(int i=0;i<level;++i) cout <<"  ";) \
       DEBUG(cout << "to memo: " << _this_rule_->name << " at: " << n << endl;) \
       memo[rule_id][n] = Memo_data(_ok_,_this_rule_,_pos_); \
       num_stored++; \
    } \
    level--; \
    return _ok_; \
}


#define RESULT() \
    if(_ok_){ \
        _this_rule_->result = [_this_rule_]()->AstNode_ptr { \
                AstNode_ptr _result_(new AstNode); \

#define RESULT1(A) \
    if(_ok_){ \
        _this_rule_->result = [_this_rule_,A]()->AstNode_ptr { \
                AstNode_ptr _result_(new AstNode); \

#define RESULT2(A,B) \
    if(_ok_){ \
        _this_rule_->result = [_this_rule_,A,B]()->AstNode_ptr { \
                AstNode_ptr _result_(new AstNode); \


/*
#define RESULT(A) \
    if(_ok_){ \
        _this_rule_->result = [_this_rule_]()->AstNode_ptr { \
                AstNode_ptr _result_(new AstNode);
*/

#define END_RESULT \
    /*DEBUG(cout << "calling result of " << _this_rule_->name << endl;) */\
    return _result_; }; }

// Rule without argument
#define LITERAL(_name, _v) \
bool _name(bool add_to_tree = true, bool do_memo = false){ \
    bool simplify = false;\
    HEADER(_name) \
    string _arg_(_v); \
    for(auto ch: _arg_){ \
        if(*_pos_==ch){ \
            _pos_++; \
        } else { \
            break; \
        } \
    } \
    if(_pos_-_old_==_arg_.size()) _ok_ = true; \
    RESULT()

#define END_LITERAL \
    END_RESULT\
    END_RULE

#define SUBRULE(n) (_this_rule_->children[n]->result())

#define NUM_SUBRULES() (_this_rule_->children.size())
/*===================================*/
/*           PREDICATES              */
/*===================================*/

// Predicate, which combines sequence of rules. If sequence fails it rewinds _pos_ to beginning
#define Comb(rule) \
    ( \
    [this,&_this_rule_]()->bool{ \
        std::string::iterator old=_pos_; \
        int old_sz = NUM_SUBRULES(); \
        bool ok = (rule);\
        if(!ok) { \
            _pos_ = old; \
            _this_rule_->children.resize(old_sz); \
        } \
        return ok; \
    }() \
    )

// Make rule optional
#define Opt(rule) (Comb(rule)||true)

// Predicate, which checks if rule matches, but doesn't advance iterator
#define Check(rule) \
    ( \
    [this,&_this_rule_]()->bool{ \
        std::string::iterator old=_pos_; \
        int old_sz = NUM_SUBRULES(); \
        bool ok = (rule);\
        _pos_ = old; \
        _this_rule_->children.resize(old_sz); \
        return ok; \
    }() \
    )

#define ZeroOrMore(rule) \
    ( \
    [this,&_this_rule_]()->bool{ \
        bool ok = true; \
        while(ok){ ok = Comb(rule); } \
        return true; \
    }() \
    )

#define OneOrMore(rule) \
    ( \
    [this,&_this_rule_]()->bool{ \
        bool ok = Comb(rule); \
        if(!ok) return false; \
        while(ok){ ok = Comb(rule); } \
        return true; \
    }() \
    )


struct Memo_data {
    Memo_data(){}
    Memo_data(bool _ok, const Parse_node_ptr& _tree, const string::iterator& _pos){
        ok = _ok;
        tree = _tree;
        pos = _pos;
    }

    bool ok; // Rule evaluation result
    Parse_node_ptr tree;  // Subtree
    string::iterator pos;  // iterator after rule completion
};



/*===================================*/
/*           THE GRAMMAR             */
/*===================================*/

class Grammar {
friend class Rule_proxy;
private:
    std::string::iterator _pos_,end,beg,last_success;
    Parse_node_ptr current_parent;
    // Rule counter
    static int rule_counter;
    // Memo table
    vector< map<int,Memo_data> > memo;

    int level; // For pretty printing
    int num_stored, num_restored;
    int max_level;

public:
    Grammar(std::string& s){
        _pos_ = beg = last_success = s.begin();
        end = s.end();
        // Initial size of memotable
        memo.reserve(100);
        level = num_stored = num_restored = max_level = 0;
    }

#ifdef _DEBUG_PARSER
    void dump(Parse_node_ptr p, int indent=0){
        for(int i=0;i<indent;++i) cout << '\t';
        cout << p->name << endl;

        for(int i=0;i<p->children.size();++i){
                dump(p->children[i],indent+1);
        }
    }
#endif

    /*===================================*/
    /*           TERMINALS               */
    /*===================================*/


    RULE(SPACE)
        while(isspace(*_pos_)) _pos_++;
        if(_old_!=_pos_) _ok_ = true;
    END_RULE

// Optional space
#define SP_() (SPACE(false,false)||true)
// Mandatory space
#define SP() (SPACE(false,false))

    RULE(UINT)
        char* e;
        int val = strtoul(&*_old_, &e, 0);
        _pos_ += e-&*_old_;
        if(_old_!=_pos_) _ok_ = true;
        SP_(); // Consume any training space if present
        RESULT1(val)
            _result_->code = TOK_UINT;
            _result_->children.push_back(val);
        END_RESULT
    END_RULE

    RULE(INT)
        char* e;
        int val = strtol(&*_old_, &e, 0);
        _pos_ += e-&*_old_;
        if(_old_!=_pos_) _ok_ = true;
        SP_(); // Consume any training space if present
        RESULT1(val)
            _result_->code = TOK_INT;
            _result_->children.push_back(val);
        END_RESULT
    END_RULE

    RULE(FLOAT)
        char* e;
        float val = strtod(&*_old_, &e);
        _pos_ += e-&*_old_;
        if(_old_!=_pos_) _ok_ = true;
        SP_(); // Consume any training space if present
        RESULT1(val)
            _result_->code = TOK_FLOAT;
            _result_->children.push_back(val);
        END_RESULT
    END_RULE


    /*===================================*/
    /*           LITERALS                */
    /*===================================*/

    LITERAL(PLUS,"+")
        _result_->code = TOK_PLUS;
    END_LITERAL

    LITERAL(MINUS,"-")
        _result_->code = TOK_MINUS;
    END_LITERAL

    LITERAL(STAR,"*")
        _result_->code = TOK_MULT;
    END_LITERAL

    LITERAL(SLASH,"/")
        _result_->code = TOK_DIV;
    END_LITERAL

    LITERAL(CAP,"^")
        _result_->code = TOK_POWER;
    END_LITERAL

    LITERAL(DOUBLE_STAR,"**")
        _result_->code = TOK_POWER;
    END_LITERAL

    LITERAL(LPAREN,"(")
    END_LITERAL

    LITERAL(RPAREN,")")
    END_LITERAL

    LITERAL(X,"x")
        _result_->code = TOK_X;
    END_LITERAL

    LITERAL(Y,"y")
        _result_->code = TOK_Y;
    END_LITERAL

    LITERAL(Z,"z")
        _result_->code = TOK_Z;
    END_LITERAL

    LITERAL(BETA,"beta")
        _result_->code = TOK_BETA;
    END_LITERAL

    LITERAL(OCCUPANCY,"occupancy")
        _result_->code = TOK_OCC;
    END_LITERAL

    LITERAL(DIST1,"distance")
    END_LITERAL

    LITERAL(DIST2,"dist")
    END_LITERAL

    LITERAL(POINT,"point")
    END_LITERAL

    LITERAL(VECTOR,"vector")
    END_LITERAL

    LITERAL(PLANE,"plane")
    END_LITERAL

    LITERAL(OR,"or")
        _result_->code = TOK_OR;
    END_LITERAL

    LITERAL(AND,"and")
        _result_->code = TOK_AND;
    END_LITERAL

    LITERAL(ALL,"all")
        _result_->code = TOK_ALL;
    END_LITERAL

    LITERAL(EQ,"==")
        _result_->code = TOK_EQ;
    END_LITERAL

    LITERAL(NEQ,"!=")
        _result_->code = TOK_NEQ;
    END_LITERAL

    LITERAL(GEQ,">=")
        _result_->code = TOK_GEQ;
    END_LITERAL

    LITERAL(LEQ,"<=")
        _result_->code = TOK_LEQ;
    END_LITERAL

    LITERAL(GT,">")
        _result_->code = TOK_GT;
    END_LITERAL

    LITERAL(LT,"<")
        _result_->code = TOK_LT;
    END_LITERAL

    LITERAL(NOT,"not")
    END_LITERAL

    LITERAL(WITHIN_,"within")
    END_LITERAL

    LITERAL(OF,"of")
    END_LITERAL

    LITERAL(PBC_ON1,"pbc")
        _result_->code = TOK_UINT;
        _result_->children.push_back(1);
    END_LITERAL

    LITERAL(PBC_ON2,"periodic")
        _result_->code = TOK_UINT;
        _result_->children.push_back(1);
    END_LITERAL

    LITERAL(PBC_OFF1,"nopbc")
        _result_->code = TOK_UINT;
        _result_->children.push_back(0);
    END_LITERAL

    LITERAL(PBC_OFF2,"noperiodic")
        _result_->code = TOK_UINT;
        _result_->children.push_back(0);
    END_LITERAL

    LITERAL(BY,"by")
    END_LITERAL

    LITERAL(TO,"to")
    END_LITERAL

    LITERAL(RESIDUE,"residue")
    END_LITERAL

    LITERAL(NAME,"name")
        _result_->code = TOK_NAME;
    END_LITERAL

    LITERAL(RESNAME,"resname")
        _result_->code = TOK_RESNAME;
    END_LITERAL

    LITERAL(TAG,"tag")
        _result_->code = TOK_TAG;
    END_LITERAL

    LITERAL(CHAIN,"chain")
        _result_->code = TOK_CHAIN;
    END_LITERAL

    LITERAL(RESID,"resid")
        _result_->code = TOK_RESID;
    END_LITERAL

    LITERAL(RESINDEX,"resindex")
        _result_->code = TOK_RESINDEX;
    END_LITERAL

    LITERAL(INDEX,"index")
        _result_->code = TOK_INDEX;
    END_LITERAL

    /*===================================*/
    /*           NON-TERMINALS           */
    /*===================================*/


    RULE_REDUCE(NUM_EXPR)

        _ok_ = NUM_TERM() && ZeroOrMore( (PLUS() || MINUS()) && SP_() && NUM_EXPR() );

        RESULT()
            _result_ = SUBRULE(0); // left
            if(NUM_SUBRULES()>1){
                for(int i=1; i<NUM_SUBRULES()-1; i+=2){
                    AstNode_ptr tmp = SUBRULE(i); // Operation
                    _result_.swap(tmp);
                    _result_->children.push_back(tmp); // left operand
                    _result_->children.push_back(SUBRULE(i+1)); // right operand
                }
            }
        END_RESULT
    END_RULE

    RULE_REDUCE(NUM_TERM)

        _ok_ = NUM_POWER() && ZeroOrMore( (STAR() || SLASH()) && SP_() && NUM_POWER() );

        RESULT()
            _result_ = SUBRULE(0); // left
            if(NUM_SUBRULES()>1){
                for(int i=1; i<NUM_SUBRULES()-1; i+=2){
                    AstNode_ptr tmp = SUBRULE(i);
                    _result_.swap(tmp);                   
                    _result_->children.push_back(tmp); // left operand
                    _result_->children.push_back(SUBRULE(i+1)); // right operand
                }
            }
        END_RESULT
    END_RULE


    RULE_REDUCE(NUM_POWER)

        _ok_ = NUM_FACTOR() && Opt( (CAP() || DOUBLE_STAR()) && SP_()&& NUM_FACTOR() );

        RESULT()
            _result_ = SUBRULE(0); // left
            if(NUM_SUBRULES()>1){
                AstNode_ptr tmp = SUBRULE(1);
                _result_.swap(tmp);                
                _result_->children.push_back(tmp); // left operand
                _result_->children.push_back(SUBRULE(2)); // right operand
            }
        END_RESULT
    END_RULE

    RULE_REDUCE(NUM_FACTOR)
        _ok_ = Comb( LPAREN(false) && SP_() && NUM_EXPR() && RPAREN(false) && SP_() )
                || FLOAT()                
                || Comb( X() && SP_() )
                || Comb( Y() && SP_() )
                || Comb( Z() && SP_() )
                || Comb( BETA() && SP_() )
                || Comb( OCCUPANCY() && SP_() )
                || UNARY_MINUS()
                || DIST_POINT()
                || DIST_VECTOR()
                || DIST_PLANE()
                ;

        RESULT()
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE    

    RULE(UNARY_MINUS)
        _ok_ = MINUS(false) && NUM_FACTOR();
        RESULT()
            _result_->code = TOK_UNARY_MINUS;
            _result_->children.push_back(SUBRULE(0));
        END_RESULT
    END_RULE

    RULE(DIST_POINT)
        _ok_ = (DIST1(false) || DIST2(false)) && SP()
                && POINT(false) && SP()
                && Opt(PBC())
                && FLOAT() && FLOAT() && FLOAT();
        RESULT()
        _result_->code = TOK_POINT;
        if(NUM_SUBRULES()==3){ // No pbc given
            _result_->children.push_back(SUBRULE(0));
            _result_->children.push_back(SUBRULE(1));
            _result_->children.push_back(SUBRULE(2));
            _result_->children.push_back(0); // default pbc
        } else {
            _result_->children.push_back(SUBRULE(1));
            _result_->children.push_back(SUBRULE(2));
            _result_->children.push_back(SUBRULE(3));
            _result_->children.push_back(SUBRULE(0)->children[0]); // pbc
        }
        END_RESULT
    END_RULE

    RULE(DIST_VECTOR)
        _ok_ = (DIST1(false) || DIST2(false)) && SP()
                && VECTOR(false) && SP()
                && Opt(PBC())
                && FLOAT() && FLOAT() && FLOAT() && FLOAT() && FLOAT() && FLOAT();
        RESULT()
        _result_->code = TOK_VECTOR;
        if(NUM_SUBRULES()==6){ // No pbc given
            _result_->children.push_back(SUBRULE(0));
            _result_->children.push_back(SUBRULE(1));
            _result_->children.push_back(SUBRULE(2));
            _result_->children.push_back(SUBRULE(3));
            _result_->children.push_back(SUBRULE(4));
            _result_->children.push_back(SUBRULE(5));
            _result_->children.push_back(0); // default pbc
        } else {
            _result_->children.push_back(SUBRULE(1));
            _result_->children.push_back(SUBRULE(2));
            _result_->children.push_back(SUBRULE(3));
            _result_->children.push_back(SUBRULE(4));
            _result_->children.push_back(SUBRULE(5));
            _result_->children.push_back(SUBRULE(6));
            _result_->children.push_back(SUBRULE(0)->children[0]); // pbc
        }
        END_RESULT
    END_RULE

    RULE(DIST_PLANE)
        _ok_ = (DIST1(false) || DIST2(false)) && SP()
                && PLANE(false) && SP()
                && Opt(PBC())
                && FLOAT() && FLOAT() && FLOAT() && FLOAT() && FLOAT() && FLOAT();
        RESULT()
        _result_->code = TOK_PLANE;
        if(NUM_SUBRULES()==6){ // No pbc given
            _result_->children.push_back(SUBRULE(0));
            _result_->children.push_back(SUBRULE(1));
            _result_->children.push_back(SUBRULE(2));
            _result_->children.push_back(SUBRULE(3));
            _result_->children.push_back(SUBRULE(4));
            _result_->children.push_back(SUBRULE(5));
            _result_->children.push_back(0); // default pbc
        } else {
            _result_->children.push_back(SUBRULE(1));
            _result_->children.push_back(SUBRULE(2));
            _result_->children.push_back(SUBRULE(3));
            _result_->children.push_back(SUBRULE(4));
            _result_->children.push_back(SUBRULE(5));
            _result_->children.push_back(SUBRULE(6));
            _result_->children.push_back(SUBRULE(0)->children[0]); // pbc
        }
        END_RESULT
    END_RULE

    RULE_REDUCE(LOGICAL_EXPR)

        _ok_ = LOGICAL_OPERAND() && ZeroOrMore( (OR() || AND()) && SP_() && LOGICAL_OPERAND() );

        RESULT()
            _result_ = SUBRULE(0); // left
            if(NUM_SUBRULES()>1){
                for(int i=1; i<NUM_SUBRULES()-1; i+=2){
                    AstNode_ptr tmp = SUBRULE(i);
                    _result_.swap(tmp);
                    _result_->children.push_back(tmp); // left operand
                    _result_->children.push_back(SUBRULE(i+1)); // right operand
                }
            }
        END_RESULT
    END_RULE

    RULE_REDUCE(LOGICAL_OPERAND)
        _ok_ = Comb( LPAREN(false) && SP_() && LOGICAL_EXPR() && RPAREN(false) && SP_() )
                ||
               Comb( !Check(NUM_EXPR(false) && !COMPARISON_OPERATOR(false)) && NUM_COMPARISON() )
                ||
               Comb( ALL() && SP_() )
                ||
               LOGICAL_NOT()
                ||
               WITHIN()
                ||
               BY_RESIDUE()
                ||
               KEYWORD_LIST_STR()
                ||
               KEYWORD_INT_STR()
                ;
        RESULT()
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(COMPARISON_OPERATOR)
        _ok_ = (EQ() || NEQ() || LEQ() || GEQ() || LT() || GT()) && SP_();
        RESULT()
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE_REDUCE(NUM_COMPARISON)
        _ok_ = NUM_EXPR();

        Comb( COMPARISON_OPERATOR() && NUM_EXPR() && COMPARISON_OPERATOR() && NUM_EXPR()) //chained
        ||
        Comb( COMPARISON_OPERATOR() && NUM_EXPR() ); // normal

        RESULT()
            if(NUM_SUBRULES()==1){ // single NUM_EXPR
                _result_ = SUBRULE(0);
            } else if(NUM_SUBRULES()==3){ // normal comparison
                _result_ = SUBRULE(1);
                _result_->children.push_back(SUBRULE(0));
                _result_->children.push_back(SUBRULE(2));
            } else { // chained comparison
                AstNode_ptr op1 = SUBRULE(1);
                op1->children.push_back(SUBRULE(0));
                op1->children.push_back(SUBRULE(2));
                AstNode_ptr op2 = SUBRULE(3);
                op2->children.push_back(SUBRULE(2));
                op2->children.push_back(SUBRULE(4));                
                _result_->code = TOK_AND;
                _result_->children.push_back(op1);
                _result_->children.push_back(op2);
            }
        END_RESULT;
    END_RULE    

    RULE(LOGICAL_NOT)
        _ok_ = NOT(false) && SP_() && LOGICAL_OPERAND();
        RESULT()
            _result_->code = TOK_NOT;
            _result_->children.push_back(SUBRULE(0));
        END_RESULT
    END_RULE

    RULE(WITHIN)
        _ok_ = WITHIN_(false) && SP_() && FLOAT() && SP_()
                && Opt(PBC()) && OF(false)
                && (SP()||Check(LPAREN(false))) && LOGICAL_OPERAND();

        RESULT()
        _result_->code = TOK_WITHIN;
        _result_->children.push_back(SUBRULE(0)->children[0]); // d
        if(NUM_SUBRULES()==2){ // no pbc given
            _result_->children.push_back(SUBRULE(1)); // operand
            _result_->children.push_back(0); // pbc
        } else { // with pbc
            _result_->children.push_back(SUBRULE(2)); // operand
            _result_->children.push_back(SUBRULE(1)->children[0]); // pbc
        }
        END_RESULT
    END_RULE

    RULE(PBC)
        _ok_ = (PBC_ON1() || PBC_ON2() || PBC_OFF1() || PBC_OFF2()) && SP_();
        RESULT()
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(BY_RESIDUE)
        _ok_ = BY(false) && SP() && RESIDUE(false) && SP_() && LOGICAL_OPERAND();
        RESULT()
        _result_->code = TOK_BY;
        _result_->children.push_back(SUBRULE(0));
        END_RESULT
    END_RULE

    RULE(KEYWORD_LIST_STR)
        _ok_ = STR_KEYWORD() && SP() && OneOrMore( STR()||REGEX() );

        RESULT()
            _result_ = SUBRULE(0);
            for(int i=1; i<NUM_SUBRULES(); ++i){
                _result_->children.push_back(SUBRULE(i));
            }
        END_RESULT
    END_RULE

    RULE(KEYWORD_INT_STR)
        _ok_ = INT_KEYWORD() && SP() && OneOrMore( RANGE()||UINT() && SP_() );

        RESULT()
            _result_ = SUBRULE(0);
            for(int i=1; i<NUM_SUBRULES(); ++i){
                _result_->children.push_back(SUBRULE(i));
            }
        END_RESULT
    END_RULE

    RULE(STR_KEYWORD)
        _ok_ = NAME() || RESNAME() || TAG() || CHAIN();
        RESULT()
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(INT_KEYWORD)
        _ok_ = RESID() || RESINDEX() || INDEX();
        RESULT()
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(STR)
        _ok_ = !Check( OR(false) || AND(false) );
        if(_ok_){
            while(isalnum(*_pos_) && _pos_!=end){
                _pos_++;
            }
        }

        string s;

        if(_pos_!=_old_){
            if(_ok_) s = string(_old_,_pos_);
            _ok_ = SP() || Check(RPAREN(false)) || Check(MINUS(false)) || (_pos_==end);
        } else {
            _ok_ = false;
        }

        RESULT1(s)
            _result_->code = TOK_STR;
            _result_->children.push_back(s);
        END_RESULT
    END_RULE

    RULE(REGEX)
        _ok_ = (*_pos_=='\'');
        if(_ok_){
            _pos_++;
            while(*_pos_!='\'' && _pos_!=end) _pos_++;
            if(_pos_!=_old_){
                _ok_ = (*_pos_=='\'');
                if(_ok_) _pos_++;
            }
        } else {
            _ok_ = (*_pos_=='"');
            if(_ok_){
                _pos_++;
                while(*_pos_!='"' && _pos_!=end) _pos_++;
                if(_pos_!=_old_){
                    _ok_ = (*_pos_=='"');
                    if(_ok_) _pos_++;
                }
            }
        }

        string::iterator b = _old_+1;
        string::iterator e = _pos_-1;        

        SP_(); // Consume any trailing space if present

        RESULT2(b,e)
            _result_->code = TOK_REGEX;
            _result_->children.push_back(string(b,e));
        END_RESULT
    END_RULE

    RULE(RANGE)
        _ok_ = UINT() && SP_() && (TO(false)||MINUS(false)) && SP_() && UINT();
        RESULT()
            _result_->code = TOK_TO;
            _result_->children.push_back(SUBRULE(0)->children[0]);
            _result_->children.push_back(SUBRULE(1)->children[0]);
        END_RESULT
    END_RULE


    RULE_REDUCE(START)
        _ok_ = SP_() && LOGICAL_EXPR();

        RESULT()
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AstNode_ptr run(){
        current_parent.reset(new Parse_node);
        //Parse_node_ptr p = current_parent;
        START();

#ifdef _DEBUG_PARSER
        cout << "Statistics:" << endl;
        cout << "Nesting level reached: " << max_level << endl;
        cout << "Size of memotable: " << memo.size() << endl;
        cout << "Rules stored to memotable: " << num_stored << endl;
        cout << "Rules restored from memotable: " << num_restored << endl;
        dump(current_parent);
        cout << "Lazily getting AST:" << endl;
#endif

        if(current_parent->children.size()>0 && _pos_== end)
            return current_parent->children[0]->result();
        else {            
            int n = std::distance(beg,last_success);
            stringstream ss;
            ss << "Syntax error in selection! Somewhere after position " << n << ":" << endl;
            ss << string(beg,end) << endl;
            for(int i=0;i<n-1;i++) ss<< "-";
            ss << "^" << endl;

            throw Pteros_error(ss.str());
        }        
    }

};

int Grammar::rule_counter = 0;




//===========================================================

#endif /* SELECTION_PARSER_H */

