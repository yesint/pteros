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

#define HEADER(_name) \
    if(_pos_==end) return false; \
    std::string::iterator _old_ = _pos_; \
    bool _ok_ = false; \
    Parse_node_ptr _this_rule_(new Parse_node); \
    _this_rule_->name = #_name; \
    Parse_node_ptr parent = current_parent; \
    current_parent = _this_rule_; \
    static int rule_used = 0; \
    static int rule_id = -1; \
    if(rule_used==0){ /*first call of this rule*/ \
        rule_id = rule_counter; \
        rule_counter++; \
    } \
    rule_used++; \
    if(rule_used==1){ if(memo.size()<rule_id+1) memo.resize(rule_id+1); } \
    int n = std::distance(beg,_pos_); \
    for(int i=0;i<level;++i) cout <<"  "; \
    cout << "::: " << _this_rule_->name << " id: " << rule_id << " at: " << n << endl; \
    /* check memotable */ \
    bool restored = false; \
    if(do_memo && memo[rule_id].count(n)==1){ \
        for(int i=0;i<level;++i) cout <<"  "; \
        cout << "from memo: " << _this_rule_->name << " at: " << n << endl; \
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


// Rule without argument
#define RULE(_name) \
bool _name(bool add_to_tree = true, bool do_memo = true){ \
    HEADER(_name) \


#define END_RULE \
    } /*if not restored*/ \
    if(_ok_){    \
        if(add_to_tree) { \
            for(int i=0;i<level;++i) cout <<"  "; \
            cout << "Adding node "<<_this_rule_->name << " to parent " << parent->name << endl; \
            parent->children.push_back(_this_rule_); \
        } \
    } else { \
        _pos_ = _old_; \
        _this_rule_->children.clear(); \
        _this_rule_->result = nullptr; \
    } \
    current_parent = parent; \
    /* Add to memotable */ \
    if(do_memo && rule_used>0 && !restored) { \
       Memo_data m; \
       m.ok = _ok_; \
       m.tree = _this_rule_; \
       m.pos = _pos_; \
       for(int i=0;i<level;++i) cout <<"  "; \
       cout << "to memo: " << _this_rule_->name << " at: " << n << endl; \
       memo[rule_id][n] = m; \
       num_stored++; \
    } \
    level--; \
    return _ok_; \
}

#define END_RESULT \
    cout << "calling result of " << _this_rule_->name << endl; \
    return _result_; }; }

#define END_LITERAL \
    END_RESULT\
    END_RULE


#define RESULT \
    if(_ok_){ \
        _this_rule_->result = [=]()->AstNode_ptr { \
                AstNode_ptr _result_(new AstNode);

#define RESULT_T AstNode_ptr

#define SUBRULE(n) (_this_rule_->children[n]->result())

// Rule without argument
#define LITERAL(_name, _v) \
bool _name(bool add_to_tree = true, bool do_memo = false){ \
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
    RESULT

/*===================================*/
/*           PREDICATES              */
/*===================================*/

// Predicate, which combines sequence of rules. If sequence fails it rewinds _pos_ to beginning
#define Comb(rule) \
    ( \
    [this,&_this_rule_]()->bool{ \
        std::string::iterator old=_pos_; \
        int old_sz = _this_rule_->children.size(); \
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
        int old_sz = _this_rule_->children.size(); \
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
    bool ok; // Rule evaluation result
    Parse_node_ptr tree;  // Subtree
    string::iterator pos;  // iterator after rule completion
};



/*===================================*/
/*           THE GRAMMAR             */
/*===================================*/
class Rule_proxy {
public:

    Rule_proxy(){
        used = 0;
    }

    bool is_cached(const string::iterator& pos){
        return memo.count(pos);
    }

    bool restore_cached(string::iterator& pos, Parse_node_ptr& parent){
        auto& m = memo[pos];
        // Add tree to parent
        parent->children.push_back(m.tree);
        // Advance iterator
        pos = m.pos;
        // Return saved result
        return m.ok;
    }

    void setup(Parse_node_ptr& current_parent){
        // Create new node and set global current_parent to it
        this_rule.reset(new Parse_node);
        saved_parent = current_parent;
        current_parent = this_rule;
    }

    void commit(Parse_node_ptr& current_parent){
        saved_parent->children.push_back(this_rule);
        current_parent = saved_parent;
    }

    void backtrack(string::iterator& pos, const string::iterator& old){
        pos = old;
        this_rule.reset();
    }

    void cache(bool ok, const string::iterator& pos){
        Memo_data m;
        m.ok = ok;
        m.tree = this_rule;
        m.pos = pos;
        memo[pos] = m;
    }

private:
    int used;
    map<string::iterator,Memo_data> memo;
    Parse_node_ptr this_rule, saved_parent;
};





class Grammar {
friend class Rule_proxy;
private:
    std::string::iterator _pos_,end,beg;
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
        _pos_ = beg = s.begin();
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
#define SP_ (SPACE(false,false)||true)
// Mandatory space
#define SP() (SPACE(false,false))

    /*
    ARG_RULE(LIT,char)
        if(*_pos_==_arg_) _ok_ = true;
        _pos_++;
        RESULT
            _result_->code = TOK_STR;
            _result_->children.push_back(_arg_);
        END_RESULT
    END_RULE

    ARG_RULE(LIT,string)
        for(auto ch: _arg_){
            if(*_pos_==ch){
                _pos_++;
            } else {
                break;
            }
        }
        if(_pos_-_old_==_arg_.size()) _ok_ = true;
        RESULT
            _result_->code = TOK_STR;
            _result_->children.push_back(_arg_);
        END_RESULT
    END_RULE
    */

    RULE(UINT)
        char* e;
        int val = strtoul(&*_old_, &e, 0);
        _pos_ += e-&*_old_;
        if(_old_!=_pos_) _ok_ = true;
        SP_; // Consume any training space if present
        RESULT
            int val = atoi(string(_old_,_pos_).c_str());
            _result_->code = TOK_UINT;
            _result_->children.push_back(val);
        END_RESULT
    END_RULE

    RULE(INT)
        char* e;
        int val = strtol(&*_old_, &e, 0);
        _pos_ += e-&*_old_;
        if(_old_!=_pos_) _ok_ = true;
        SP_; // Consume any training space if present
        RESULT
            int val = atoi(string(_old_,_pos_).c_str());
            _result_->code = TOK_INT;
            _result_->children.push_back(val);
        END_RESULT
    END_RULE

    RULE(FLOAT)
        char* e;
        float val = strtod(&*_old_, &e);
        _pos_ += e-&*_old_;
        if(_old_!=_pos_) _ok_ = true;
        SP_; // Consume any training space if present
        RESULT
            _result_->code = TOK_FLOAT;
            _result_->children.push_back(val);
        END_RESULT
    END_RULE

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


    RULE(NUM_EXPR)

        //_ok_ = NUM_TERM() && ZeroOrMore( (LIT('+') || LIT('-')) && SP_ && NUM_EXPR() );
        _ok_ = NUM_TERM() && ZeroOrMore( (PLUS() || MINUS()) && SP_ && NUM_EXPR() );

        RESULT
            _result_ = SUBRULE(0); // left
            if(_this_rule_->children.size()>1){
                for(int i=1; i<_this_rule_->children.size()-1; i+=2){
                    AstNode_ptr tmp = SUBRULE(i); // Operation
                    _result_.swap(tmp);
                    _result_->children.push_back(tmp); // left operand
                    _result_->children.push_back(SUBRULE(i+1)); // right operand
                }
            }
        END_RESULT
    END_RULE

    RULE(NUM_TERM)

        _ok_ = NUM_POWER() && ZeroOrMore( (STAR() || SLASH()) && SP_ && NUM_POWER() );

        RESULT
            _result_ = SUBRULE(0); // left
            if(_this_rule_->children.size()>1){
                for(int i=1; i<_this_rule_->children.size()-1; i+=2){
                    AstNode_ptr tmp = SUBRULE(i);
                    _result_.swap(tmp);                   
                    _result_->children.push_back(tmp); // left operand
                    _result_->children.push_back(SUBRULE(i+1)); // right operand
                }
            }
        END_RESULT
    END_RULE


    RULE(NUM_POWER)

        _ok_ = NUM_FACTOR() && Opt( (CAP() || DOUBLE_STAR()) && SP_&& NUM_FACTOR() );

        RESULT
            _result_ = SUBRULE(0); // left
            if(_this_rule_->children.size()>1){
                AstNode_ptr tmp = SUBRULE(1);
                _result_.swap(tmp);                
                _result_->children.push_back(tmp); // left operand
                _result_->children.push_back(SUBRULE(2)); // right operand
            }
        END_RESULT
    END_RULE

    RULE(NUM_FACTOR)
        _ok_ = Comb( LPAREN(false) && SP_ && NUM_EXPR() && RPAREN(false) && SP_ )
                || FLOAT()                
                || Comb( X() && SP_ )
                || Comb( Y() && SP_ )
                || Comb( Z() && SP_ )
                || Comb( BETA() && SP_ )
                || Comb( OCCUPANCY() && SP_ )
                || UNARY_MINUS()
                || DIST_POINT()
                || DIST_VECTOR()
                || DIST_PLANE()
                ;

        RESULT
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE    

    RULE(UNARY_MINUS)
        _ok_ = MINUS(false) && NUM_FACTOR();
        RESULT            
            _result_->code = TOK_UNARY_MINUS;
            _result_->children.push_back(SUBRULE(0));
        END_RESULT
    END_RULE

    RULE(DIST_POINT)
        _ok_ = (DIST1(false) || DIST2(false)) && SP()
                && POINT(false) && SP()
                && Opt(PBC())
                && FLOAT() && FLOAT() && FLOAT();
        RESULT
        _result_->code = TOK_POINT;
        if(_this_rule_->children.size()==3){ // No pbc given
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
        RESULT
        _result_->code = TOK_VECTOR;
        if(_this_rule_->children.size()==6){ // No pbc given
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
        RESULT
        _result_->code = TOK_PLANE;
        if(_this_rule_->children.size()==6){ // No pbc given
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

    RULE(LOGICAL_EXPR)

        _ok_ = LOGICAL_OPERAND() && ZeroOrMore( (OR() || AND()) && SP_ && LOGICAL_OPERAND() );

        RESULT
            _result_ = SUBRULE(0); // left
            if(_this_rule_->children.size()>1){
                for(int i=1; i<_this_rule_->children.size()-1; i+=2){
                    AstNode_ptr tmp = SUBRULE(i);
                    _result_.swap(tmp);
                    _result_->children.push_back(tmp); // left operand
                    _result_->children.push_back(SUBRULE(i+1)); // right operand
                }
            }
        END_RESULT
    END_RULE

    RULE(LOGICAL_OPERAND)
        _ok_ = Comb( LPAREN(false) && SP_ && LOGICAL_EXPR() && RPAREN(false) && SP_ )
                ||
               Comb( !Check(NUM_EXPR(false) && !COMPARISON_OPERATOR(false)) && NUM_COMPARISON() )
                ||
               Comb( ALL() && SP_ )
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
        RESULT
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(COMPARISON_OPERATOR)
        _ok_ = (EQ() || NEQ() || LEQ() || GEQ() || LT() || GT()) && SP_;
        RESULT                        
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(NUM_COMPARISON)
        _ok_ = NUM_EXPR();

        Comb( COMPARISON_OPERATOR() && NUM_EXPR() && COMPARISON_OPERATOR() && NUM_EXPR()) //chained
        ||
        Comb( COMPARISON_OPERATOR() && NUM_EXPR() ); // normal

        RESULT
            if(_this_rule_->children.size()==1){ // single NUM_EXPR
                _result_ = SUBRULE(0);
            } else if(_this_rule_->children.size()==3){ // normal comparison
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
        _ok_ = NOT(false) && SP_ && LOGICAL_OPERAND();
        RESULT            
            _result_->code = TOK_NOT;
            _result_->children.push_back(SUBRULE(0));
        END_RESULT
    END_RULE

    RULE(WITHIN)
        _ok_ = WITHIN_(false) && SP_ && FLOAT() && SP_
                && Opt(PBC()) && OF(false)
                && (SP()||Check(LPAREN(false))) && LOGICAL_OPERAND();

        RESULT        
        _result_->code = TOK_WITHIN;
        _result_->children.push_back(SUBRULE(0)->children[0]); // d
        if(_this_rule_->children.size()==2){ // no pbc given
            _result_->children.push_back(SUBRULE(1)); // operand
            _result_->children.push_back(0); // pbc
        } else { // with pbc
            _result_->children.push_back(SUBRULE(2)); // operand
            _result_->children.push_back(SUBRULE(1)->children[0]); // pbc
        }
        END_RESULT
    END_RULE

    RULE(PBC)
        _ok_ = (PBC_ON1() || PBC_ON2() || PBC_OFF1() || PBC_OFF2()) && SP_;
        RESULT            
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(BY_RESIDUE)
        _ok_ = BY(false) && SP() && RESIDUE(false) && SP_ && LOGICAL_OPERAND();
        RESULT        
        _result_->code = TOK_BY;
        _result_->children.push_back(SUBRULE(0));
        END_RESULT
    END_RULE

    RULE(KEYWORD_LIST_STR)
        _ok_ = STR_KEYWORD() && SP() && OneOrMore( STR()||REGEX() );

        RESULT
            _result_ = SUBRULE(0);
            for(int i=1; i<_this_rule_->children.size(); ++i){
                _result_->children.push_back(SUBRULE(i));
            }
        END_RESULT
    END_RULE

    RULE(KEYWORD_INT_STR)
        _ok_ = INT_KEYWORD() && SP() && OneOrMore( RANGE()||UINT() && SP_ );

        RESULT
            _result_ = SUBRULE(0);
            for(int i=1; i<_this_rule_->children.size(); ++i){
                _result_->children.push_back(SUBRULE(i));
            }
        END_RESULT
    END_RULE

    RULE(STR_KEYWORD)
        _ok_ = NAME() || RESNAME() || TAG() || CHAIN();
        RESULT            
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(INT_KEYWORD)
        _ok_ = RESID() || RESINDEX() || INDEX();
        RESULT            
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

        RESULT            
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

        SP_; // Consume any trailing space if present

        RESULT            
            _result_->code = TOK_REGEX;
            _result_->children.push_back(string(b,e));
        END_RESULT
    END_RULE

    RULE(RANGE)
        _ok_ = UINT() && SP_ && (TO(false)||MINUS(false)) && SP_ && UINT();
        RESULT            
            _result_->code = TOK_TO;
            _result_->children.push_back(SUBRULE(0)->children[0]);
            _result_->children.push_back(SUBRULE(1)->children[0]);
        END_RESULT
    END_RULE


    RULE(START)
        _ok_ = SP_ && LOGICAL_EXPR();

        RESULT
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AstNode_ptr run(){
        current_parent.reset(new Parse_node);
        Parse_node_ptr p = current_parent;
        START();

        cout << "Statistics:" << endl;
        cout << "Nesting level reached: " << max_level << endl;
        cout << "Size of memotable: " << memo.size() << endl;
        cout << "Rules stored to memotable: " << num_stored << endl;
        cout << "Rules restored from memotable: " << num_restored << endl;

#ifdef _DEBUG_PARSER
        dump(p);
        cout << "Lazily getting AST:" << endl;
#endif
        if(p->children.size()>0 && _pos_== end)
            return p->children[0]->result();
        else {
            cout << "Syntax error at " << *_pos_ << endl;
            return nullptr;
        }        
    }


    bool rule_a(){
        static Rule_proxy* proxy = nullptr;
        if(!proxy){
            rules.push_back(Rule_proxy());
            proxy = &rules.back();
        }
        if(proxy->is_cached(_pos_)){
            return proxy->restore_cached(_pos_,current_parent);
        }

        // If we are here than not chached
        proxy->setup(current_parent);

        // Variables to be used in body
        bool _ok_ = false;
        string::iterator _old_ = _pos_;

        //////////////////
        // Here is a body
        //////////////////

        // Footer
        if(_ok_){
            proxy->commit(current_parent);
        } else {
            proxy->backtrack(_pos_,_old_);
        }

        // Add to memotable if needed
        proxy->cache(_ok_,_pos_);
        return _ok_;
    }

private:
    vector<Rule_proxy> rules;
};

int Grammar::rule_counter = 0;




//===========================================================


#endif /* SELECTION_PARSER_H */
