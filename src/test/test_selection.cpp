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

#include <string>
#include "pteros/analysis/options.h"
#include <Eigen/Core>
#include "pteros/core/pteros_error.h"
#include "pteros/core/selection.h"

#include "pteros/core/peg_parser.h"

#include <boost/variant.hpp>
#include <functional>
#include <list>

using namespace std;
using namespace pteros;
using namespace Eigen;
/*
struct Parser {
    void parse(){
        input_data INP;
        //INP.string_to_parse = "-( 2.3e-8+ (x - 3.014 ) *1.2)-(3+4- -4*5.5-(2+-3.1)) ";
        INP.string_to_parse = "(3.1+1)**5.0+2*2<x";
        INP.cur_buf_pos = 0;
        GREG g;
        yyinit(&g);
        g.data = &INP;
        while (yyparse(&g));
        if(g.begin != g.end) throw Pteros_error("Syntax error at pos "+to_string(g.begin));

        yydeinit(&g);
        INP.root->dump();
    }
};
*/

//===========================================================

typedef std::function<AstNode_ptr()> result_t;

class Parse_node {
public:
    std::string name; // Rule name
    result_t result; // result of the rule (code returning AstNode_ptr)
    std::vector<std::shared_ptr<Parse_node>> children; // Sub-rules
};

typedef std::shared_ptr<Parse_node> Parse_node_ptr;

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

// Modifier, which sais not to add this rule to the tree
#define SKIP add_to_tree = false;

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

#define RESULT \
    if(_ok_){ \
        _this_rule_->result = [=]()->AstNode_ptr { \
                AstNode_ptr _result_(new AstNode);

#define END_RESULT \
    cout << "calling result of " << _this_rule_->name << " " << _arg_ << endl; \
    return _result_; }; }

#define RESULT_T AstNode_ptr

#define SUBRULE(n) (_this_rule_->children[n]->result())

// Make rule optional
#define opt(rule) (rule||true)

// Predicate, which checks if rule matches, but doesn't advance iterator
#define check(rule) \
    ( \
    [this]()->bool{ \
        string::iterator old=_pos_; \
        bool ok = rule;\
        _pos_ = old; \
        return ok; \
    }() \
    )

class Grammar {
public:
    Grammar(std::string s){
        str = s;
        _pos_ = str.begin();
        end = str.end();
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
        RESULT            
            _result_->code = TOK_FLOAT;
            _result_->children.push_back(val);
        END_RESULT
    END_RULE

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
                || ( FLOAT() && opt(SP()) )
                || ( INT() && opt(SP()) )
                || ( X() && opt(SP()) )
                || ( Y() && opt(SP()) )
                || ( Z() && opt(SP()) )
                || ( BETA() && opt(SP()) )
                || ( OCCUPANCY() && opt(SP()) )
                || UNARY_MINUS()
                ;

        RESULT
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    RULE(X)
        _ok_ = LIT('x',false);
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_X;
        END_RESULT
    END_RULE

    RULE(Y)
        _ok_ = LIT('y',false);
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_Y;
        END_RESULT
    END_RULE

    RULE(Z)
        _ok_ = LIT('z',false);
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_Z;
        END_RESULT
    END_RULE

    RULE(BETA)
        _ok_ = LIT("beta",false);
        RESULT
            _result_.reset(new AstNode);
            _result_->code = TOK_BETA;
        END_RESULT
    END_RULE

    RULE(OCCUPANCY)
        _ok_ = LIT("occupancy",false);
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

    RULE(START)
        _ok_ = opt(SP()) && NUM_EXPR();

        RESULT
            _result_ = SUBRULE(0);
        END_RESULT
    END_RULE

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AstNode_ptr run(){
        current_parent.reset(new Parse_node);
        Parse_node_ptr p = current_parent;
        START();
        cout << "Lazily getting AST:" << endl;
        if(p->children.size()>0)
            return p->children[0]->result();
        else {
            cout << "Syntax error at " << *_pos_ << endl;
            return nullptr;
        }
    }

private:
    std::string::iterator _pos_,end;
    string str;
    Parse_node_ptr current_parent;
};

//===========================================================

int main(int argc, char** argv)
{

    try{        

        Grammar g(" -x* -(y/x^2.2)");
        AstNode_ptr res = g.run();
        if(res) res->dump();

        //cout << (boost::get<Parse_tree_ptr>(p->children.front())) << endl;

        //std::shared_ptr<Parser> p(new Parser);
        //p->parse();

//        System s("/home/semen/work/Projects/asymmetric_bilayer/for-diamonds/hexagonal/2x2.gro");
  //      Selection sel(s,"(y+4)<(x+2) or (1-z)>x");


        /*
        string str("--trajectory["
                   "initial_structure.pdb "
                   "traj.xtc r1/traj.xtc r1/traj.part0002.xtc r1/traj.part0003.xtc "
                   "--first_frame 0"
                   "--last_frame 100 "
                   "--log_interval 2 "
                   "] "
                   "--task rms[--selection \"name CA\" --unwrap 0.2] "
                   "--task rms[--selection \"protein\" --unwrap 0.4] "
                   "--dump_input dump");
        cout << str << endl;

        Options_tree opt;
        opt.from_command_line(str);
        for(auto o: opt.get_options("task")){
            cout << o->get_value<string>("selection") << endl;
            cout << o->get_value<double>("unwrap") << endl;
        }

        System s;


        Options toplevel;
        vector<Options> tasks;

        parse_command_line(argc,argv,toplevel,"task",tasks);
        //parse_command_line(argc,argv,toplevel);
        toplevel.debug();
        cout << "------" << endl;
        for(auto& t: tasks){
            t.debug();
            cout << "------" << endl;
        }

        vector<float> v = toplevel("tramvay","3.14159 42 -4.5").as_floats();
        v = toplevel("tramvay").as_floats();
        for(auto a: v) cout << a << endl;

*/
    } catch(const Pteros_error& e){ e.print(); }

}

