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

// Rule with argument
#define ARG_RULE(_name, _arg_t) \
bool _name(_arg_t _arg_, bool add_to_tree = true){ \
    if(_pos_==end) return false; \
    std::string::iterator _old_ = _pos_; \
    bool _ok_ = false; \
    Parse_node_ptr _this_rule_(new Parse_node); \
    _this_rule_->name = #_name; \
    Parse_node_ptr parent = current_parent; \
    current_parent = _this_rule_; \

// Rule without argument
#define RULE(_name) \
bool _name(bool add_to_tree = true){ \
    if(_pos_==end) return false; \
    bool _arg_ = false; /* dummy _arg_ */\
    std::string::iterator _old_ = _pos_; \
    bool _ok_ = false; \
    Parse_node_ptr _this_rule_(new Parse_node); \
    _this_rule_->name = #_name; \
    Parse_node_ptr parent = current_parent; \
    current_parent = _this_rule_; \


#define END_RULE \
    if(_ok_){    \
        if(add_to_tree) { \
            cout << "Adding node "<<_this_rule_->name << " " << _arg_ << " to parent " << parent->name << endl; \
            parent->children.push_back(_this_rule_); \
        } \
    } else { \
        _pos_ = _old_; \
        _this_rule_->children.clear(); \
        _this_rule_->result = nullptr; \
    } \
    current_parent = parent; \
    return _ok_; \
}

#define END_RESULT \
    cout << "calling result of " << _this_rule_->name << " " << _arg_ << endl; \
    return _result_; }; }


#else //DEBUG

// Rule with argument
#define ARG_RULE(_name, _arg_t) \
bool _name(_arg_t _arg_, bool add_to_tree = true){ \
    if(_pos_==end) return false; \
    std::string::iterator _old_ = _pos_; \
    bool _ok_ = false; \
    Parse_node_ptr _this_rule_(new Parse_node); \
    Parse_node_ptr parent = current_parent; \
    current_parent = _this_rule_;

// Rule without argument
#define RULE(_name) \
bool _name(bool add_to_tree = true){ \
    if(_pos_==end) return false; \
    bool _arg_ = false; /* dummy _arg_ */\
    std::string::iterator _old_ = _pos_; \
    bool _ok_ = false; \
    Parse_node_ptr _this_rule_(new Parse_node); \
    Parse_node_ptr parent = current_parent; \
    current_parent = _this_rule_;

#define END_RULE \
    if(_ok_){    \
        if(add_to_tree) { \
            parent->children.push_back(_this_rule_); \
        } \
    } else { \
        _pos_ = _old_; \
        _this_rule_->children.clear(); \
        _this_rule_->result = nullptr; \
    } \
    current_parent = parent; \
    return _ok_; \
}

#define END_RESULT \
    return _result_; }; }


#endif //DEBUG


// Modifier, which sais not to add this rule to the tree
#define SKIP add_to_tree = false;

#define RESULT \
    if(_ok_){ \
        _this_rule_->result = [=]()->AstNode_ptr { \
                AstNode_ptr _result_(new AstNode);

#define RESULT_T AstNode_ptr

#define SUBRULE(n) (_this_rule_->children[n]->result())

/*===================================*/
/*           PREDICATES              */
/*===================================*/

// Make rule optional
#define opt(rule) (rule||true)

// Predicate, which checks if rule matches, but doesn't advance iterator
#define check(rule) \
    ( \
    [this]()->bool{ \
        std::string::iterator old=_pos_; \
        bool ok = rule;\
        _pos_ = old; \
        return ok; \
    }() \
    )

class Grammar {

private:
    std::string::iterator _pos_,end;
    string str;
    Parse_node_ptr current_parent;

public:
    Grammar(std::string s){
        str = s;
        _pos_ = str.begin();
        end = str.end();
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


    RULE(SP) SKIP
        while(isspace(*_pos_)) _pos_++;
        if(_old_!=_pos_) _ok_ = true;
    END_RULE

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

    RULE(UINT)
        char* e;
        int val = strtoul(&*_old_, &e, 0);
        _pos_ += e-&*_old_;
        if(_old_!=_pos_) _ok_ = true;
        opt(SP()); // Consume any training space if present
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
        opt(SP()); // Consume any training space if present
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
        opt(SP()); // Consume any training space if present
        RESULT
            _result_->code = TOK_FLOAT;
            _result_->children.push_back(val);
        END_RESULT
    END_RULE

    /*===================================*/
    /*           NON-TERMINALS           */
    /*===================================*/


    RULE(NUM_EXPR)
        _ok_ = NUM_TERM();
        bool ok = true;
        while(ok){
             ok = (LIT('+') || LIT('-')) && opt(SP()) && NUM_EXPR();
        }

        RESULT
            _result_ = SUBRULE(0); // left
            if(_this_rule_->children.size()>1){
                for(int i=1; i<_this_rule_->children.size()-1; i+=2){
                    auto tmp = AstNode_ptr(new AstNode);
                    _result_.swap(tmp);

                    if(boost::get<char>(SUBRULE(i)->children[0])=='+')
                        _result_->code = TOK_PLUS;
                    else
                        _result_->code = TOK_MINUS;

                    _result_->children.push_back(tmp); // left operand
                    _result_->children.push_back(SUBRULE(i+1)); // right operand
                }
            }
        END_RESULT
    END_RULE

    RULE(NUM_TERM)
        _ok_ = NUM_POWER();
        bool ok = true;
        while(ok){
             ok = (LIT('*') || LIT('/')) && opt(SP()) && NUM_POWER();
        }

        RESULT
            _result_ = SUBRULE(0); // left
            if(_this_rule_->children.size()>1){
                for(int i=1; i<_this_rule_->children.size()-1; i+=2){
                    auto tmp = AstNode_ptr(new AstNode);
                    _result_.swap(tmp);

                    if(boost::get<char>(SUBRULE(i)->children[0])=='*')
                        _result_->code = TOK_MULT;
                    else
                        _result_->code = TOK_DIV;

                    _result_->children.push_back(tmp); // left operand
                    _result_->children.push_back(SUBRULE(i+1)); // right operand
                }
            }
        END_RESULT
    END_RULE


    RULE(NUM_POWER)
        _ok_ = NUM_FACTOR();
        if(_ok_){
              (LIT("^",false) || LIT("**",false)) && opt(SP())&& NUM_FACTOR();
        }

        RESULT
            _result_ = SUBRULE(0); // left
            if(_this_rule_->children.size()>1){
                auto tmp = AstNode_ptr(new AstNode);
                _result_.swap(tmp);
                _result_->code = TOK_POWER;
                _result_->children.push_back(tmp); // left operand
                _result_->children.push_back(SUBRULE(1)); // right operand
            }
        END_RESULT
    END_RULE

    RULE(NUM_FACTOR)
        _ok_ = LIT('(',false) && opt(SP()) && NUM_EXPR() && LIT(')',false) && opt(SP())
                || FLOAT()
                //|| ( INT() && opt(SP()) )
                || X()
                || Y()
                || Z()
                || BETA()
                || OCCUPANCY()
                || UNARY_MINUS()
                || DIST_POINT()
                || DIST_VECTOR()
                || DIST_PLANE()
                ;

        RESULT
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(X)
        _ok_ = LIT('x',false);
        opt(SP()); // Consume any training space if present
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_X;
        END_RESULT
    END_RULE

    RULE(Y)
        _ok_ = LIT('y',false);
        opt(SP()); // Consume any training space if present
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_Y;
        END_RESULT
    END_RULE

    RULE(Z)
        _ok_ = LIT('z',false);
        opt(SP()); // Consume any training space if present
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_Z;
        END_RESULT
    END_RULE

    RULE(BETA)
        _ok_ = LIT("beta",false);
        opt(SP()); // Consume any training space if present
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_BETA;
        END_RESULT
    END_RULE

    RULE(OCCUPANCY)
        _ok_ = LIT("occupancy",false);
        opt(SP()); // Consume any training space if present
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_OCC;
        END_RESULT
    END_RULE

    RULE(UNARY_MINUS)
        _ok_ = LIT('-',false) && NUM_FACTOR();
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_UNARY_MINUS;
            _result_->children.push_back(SUBRULE(0));
        END_RESULT
    END_RULE

    RULE(DIST_POINT)
        _ok_ = (LIT("distance",false) || LIT("dist",false)) && SP()
                && LIT("point",false) && SP()
                && opt(PBC())
                && FLOAT() && FLOAT() && FLOAT();
        RESULT
        _result_.reset(new AstNode);
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
        _ok_ = (LIT("distance",false) || LIT("dist",false)) && SP()
                && LIT("vector",false) && SP()
                && opt(PBC())
                && FLOAT() && FLOAT() && FLOAT() && FLOAT() && FLOAT() && FLOAT();
        RESULT
        _result_.reset(new AstNode);
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
        _ok_ = (LIT("distance",false) || LIT("dist",false)) && SP()
                && LIT("plane",false) && SP()
                && opt(PBC())
                && FLOAT() && FLOAT() && FLOAT() && FLOAT() && FLOAT() && FLOAT();
        RESULT
        _result_.reset(new AstNode);
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
        _ok_ = LOGICAL_OPERAND();
        bool ok = true;
        while(ok){
             ok = (LIT("or") || LIT("and")) && opt(SP()) && LOGICAL_OPERAND();
        }

        RESULT
            _result_ = SUBRULE(0); // left
            if(_this_rule_->children.size()>1){
                for(int i=1; i<_this_rule_->children.size()-1; i+=2){
                    auto tmp = AstNode_ptr(new AstNode);
                    _result_.swap(tmp);

                    if(boost::get<string>(SUBRULE(i)->children[0])=="or")
                        _result_->code = TOK_OR;
                    else
                        _result_->code = TOK_AND;

                    _result_->children.push_back(tmp); // left operand
                    _result_->children.push_back(SUBRULE(i+1)); // right operand
                }
            }
        END_RESULT
    END_RULE

    RULE(LOGICAL_OPERAND)
        _ok_ = ( LIT('(',false) && opt(SP()) && LOGICAL_EXPR() && LIT(')',false) && opt(SP()) )
                ||
               ( !check(NUM_EXPR(false) && !COMPARISON_OPERATOR(false)) && NUM_COMPARISON() )
                ||
               ALL()
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
        _ok_ = (LIT("==") || LIT("!=") || LIT("<") || LIT(">") || LIT("<=") || LIT(">=")) && opt(SP());
        RESULT
            _result_.reset(new AstNode);
            string s = boost::get<string>(SUBRULE(0)->children[0]);
            if     (s=="==") _result_->code = TOK_EQ;
            else if(s=="!=") _result_->code = TOK_NEQ;
            else if(s=="<")  _result_->code = TOK_LT;
            else if(s==">")  _result_->code = TOK_GT;
            else if(s=="<=") _result_->code = TOK_LEQ;
            else if(s==">=") _result_->code = TOK_GEQ;
        END_RESULT
    END_RULE

    RULE(NUM_COMPARISON)
        _ok_ = NUM_EXPR();
        ( COMPARISON_OPERATOR() && NUM_EXPR() && COMPARISON_OPERATOR() && NUM_EXPR()) //chained
        ||
        ( COMPARISON_OPERATOR() && NUM_EXPR() ); // normal
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
                _result_.reset(new AstNode);
                _result_->code = TOK_AND;
                _result_->children.push_back(op1);
                _result_->children.push_back(op2);
            }
        END_RESULT;
    END_RULE

    RULE(ALL)
        _ok_ = LIT("all",false) && opt(SP());
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_ALL;
        END_RESULT
    END_RULE

    RULE(LOGICAL_NOT)
        _ok_ = LIT("not",false) && opt(SP()) && LOGICAL_OPERAND();
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_NOT;
            _result_->children.push_back(SUBRULE(0));
        END_RESULT
    END_RULE

    RULE(WITHIN)
        _ok_ = LIT("within",false) && opt(SP()) && FLOAT() && opt(SP())
                && opt(PBC()) && LIT("of",false)
                && (SP()||check(LIT('(',false))) && LOGICAL_OPERAND();

        RESULT
        _result_.reset(new AstNode);
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
        _ok_ = (LIT("pbc") || LIT("periodic") || LIT("nopbc") || LIT("noperiodic")) && opt(SP());
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_INT;
            string s = boost::get<string>(SUBRULE(0)->children[0]);
            if(s=="pbc" || s=="periodic"){
                _result_->children.push_back(1);
            } else {
                _result_->children.push_back(0);
            }
        END_RESULT
    END_RULE

    RULE(BY_RESIDUE)
        _ok_ = LIT("by",false) && SP() && LIT("residue",false) && opt(SP()) && LOGICAL_OPERAND();
        RESULT
        _result_.reset(new AstNode);
        _result_->code = TOK_BY;
        _result_->children.push_back(SUBRULE(0));
        END_RESULT
    END_RULE

    RULE(KEYWORD_LIST_STR)
        _ok_ = STR_KEYWORD() && SP() && (STR()||REGEX());
        while( STR()||REGEX() );
        RESULT
            _result_ = SUBRULE(0);
            for(int i=1; i<_this_rule_->children.size(); ++i){
                _result_->children.push_back(SUBRULE(i));
            }
        END_RESULT
    END_RULE

    RULE(KEYWORD_INT_STR)
        _ok_ = INT_KEYWORD() && SP() && (RANGE()||UINT() && opt(SP()));
        while( (RANGE()||UINT()) && opt(SP()) );
        RESULT
            _result_ = SUBRULE(0);
            for(int i=1; i<_this_rule_->children.size(); ++i){
                _result_->children.push_back(SUBRULE(i));
            }
        END_RESULT
    END_RULE

    RULE(STR_KEYWORD)
        _ok_ = LIT("name") || LIT("resname") || LIT("tag") || LIT("chain");
        RESULT
            _result_.reset(new AstNode);
            string s = boost::get<string>(SUBRULE(0)->children[0]);
            if     (s=="name") _result_->code = TOK_NAME;
            else if(s=="resname") _result_->code = TOK_RESNAME;
            else if(s=="tag") _result_->code = TOK_TAG;
            else if(s=="chain") _result_->code = TOK_CHAIN;
        END_RESULT
    END_RULE

    RULE(INT_KEYWORD)
        _ok_ = LIT("resid") || LIT("resindex") || LIT("index");
        RESULT
            _result_.reset(new AstNode);
            string s = boost::get<string>(SUBRULE(0)->children[0]);
            if     (s=="resid") _result_->code = TOK_RESID;
            else if(s=="resindex") _result_->code = TOK_RESINDEX;
            else if(s=="index") _result_->code = TOK_INDEX;
        END_RESULT
    END_RULE

    RULE(STR)
        _ok_ = !check(LIT("or",false)||LIT("and",false));
        if(_ok_){
            while(isalnum(*_pos_) && _pos_!=end){
                _pos_++;
            }
        }
        if(_pos_!=_old_){
            _ok_ = SP() || check(LIT(')',false)) || check(LIT('-',false)) || (_pos_==end);
        } else {
            _ok_ = false;
        }

        string s;
        if(_ok_) s = string(_old_,_pos_);

        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_STR;
            _result_->children.push_back(s);
        END_RESULT
    END_RULE

    RULE(REGEX)
        _ok_ = LIT('\'',false);
        if(_ok_){
            while(*_pos_!='\'' && _pos_!=end) _pos_++;
            if(_pos_!=_old_) _ok_ = LIT('\'',false);
        } else {
            _ok_ = LIT('"',false);
            if(_ok_){
                while(*_pos_!='"' && _pos_!=end) _pos_++;
                if(_pos_!=_old_) _ok_ = LIT('"',false);
            }
        }

        string::iterator b = _old_+1;
        string::iterator e = _pos_-1;

        opt(SP()); // Consume any trailing space if present

        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_REGEX;
            _result_->children.push_back(string(b,e));
        END_RESULT
    END_RULE

    RULE(RANGE)
        _ok_ = UINT() && opt(SP()) && (LIT("to",false)||LIT('-',false)) && opt(SP()) && UINT();
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_TO;
            _result_->children.push_back(SUBRULE(0)->children[0]);
            _result_->children.push_back(SUBRULE(1)->children[0]);
        END_RESULT
    END_RULE


    RULE(START)
        _ok_ = opt(SP()) && LOGICAL_EXPR();

        RESULT
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AstNode_ptr run(){
        current_parent.reset(new Parse_node);
        Parse_node_ptr p = current_parent;
        START();

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


};

//===========================================================


#endif /* SELECTION_PARSER_H */
