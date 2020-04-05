//auto INTEGER = g.rule() << lit('-')('?') << any_char("0123456789")('+');
//auto FLOAT = g.rule() << INTEGER << (any_char("eE") << INTEGER)('?')

#include <list>
#include <map>
#include <string>
#include <iostream>
#include <functional>
#include <memory>
#include <variant>


using namespace std;

class Grammar;

struct Ast_node {
    int id;
    list<shared_ptr<Ast_node>> nodes;
    std::variant<string,int,float> value;
    Ast_node(int i): id(i) {}

    string to_str(Grammar* gr, string prefix="");
};


// Classes to match different rules
struct Rule_result {
    bool ok;
    string::const_iterator first,last;
    // AST
    shared_ptr<Ast_node> ast;

    Rule_result(): ok(false), first(nullptr), last(nullptr), ast(nullptr) {}
};




class Rule_holder;

// Rule class
class Rule {
    friend class Grammar;
public:
    enum class Kind {none,sequence,alternative,literal};

    Rule(): kind(Kind::none), name(""), id(-1), enable_ast(false) {}
    Rule(string s): kind(Kind::none), name(s), id(-1), enable_ast(false) {}

    string name;
    int id; // Rule tag assigned to AST nodes
    Kind kind;
    Grammar* grammar;
    list<Rule_holder> rules;
    bool enable_ast;

    // Set grammar recursively in all children
    void set_grammar(Grammar* gr);

    string to_str(string prefix="");

    virtual Rule_result do_match(bool gen_ast) = 0;
    Rule_result match(bool gen_ast);
};

// Structure which holds a polymorphic pointer to rule
// and provides operators for PEG operators
class Rule_holder {
public:
    Rule_holder(Rule* r): rule(r) {}
    Rule* operator->(){return rule;}
    Rule* get(){return rule;}

    Rule_holder operator<=(const Rule_holder& rhs){
        rule->rules.push_back(rhs);
        return *this;
    }

    Rule_holder operator<<(const Rule_holder& rhs);

    Rule_holder operator|(const Rule_holder& rhs);

    Rule_holder operator[](char c);

    Rule_holder operator!();


private:
    Rule* rule;
};


//---------------------------------------------------------
// Grammar class

#define RULE(name) \
    const int NODE_##name = __COUNTER__; \
    Rule_holder name = add_rule(new Rule_id(this,#name,NODE_##name))

class Grammar {
public:

    Grammar(): skip_whitespaces(true) {}

    void set_input(const string& inp){
        input = inp;
        cur = input.begin();
        end_input = input.end();
    }

    Rule_result parse(){
        // If asked to skip whitespace, skip in the beginning
        while(isspace(*cur)) cur++;

        top_rule->set_grammar(this);
        auto res = top_rule->match(true);
        if(res.ok){
            // Check if all input consumed
            if(cur != end_input) res.ok = false;
        }
        if(!res.ok){
            cout << "Parsing failed here:" << endl;
            cout << input << endl;
            cout << string(res.last-input.begin()-1,'_') << "^" << endl;
        }
        return res;
    }

    string to_string(){
        return top_rule->to_str();
    }

    // Input string
    string input;
    // Current position
    string::const_iterator cur;
    // End of input
    string::const_iterator end_input;

    bool skip_whitespaces;

    Rule* top_rule;

    map<int,string> id_names;

    // List of rules
    list<unique_ptr<Rule>> rules;
    // Helper method to create new rule
    // Called as grammar->add_rule(new Rule_derived())
    // Returns a holder for newly added rule
    Rule_holder add_rule(Rule* ptr){
        ptr->grammar = this; // important to set grammar for all rules!
        rules.emplace_back(ptr);
        return Rule_holder(rules.rbegin()->get());
    }

    // Helper method to replace the rule

    // Predefined rules
    Rule_holder lit(const string& s);
    Rule_holder ch(const string& s);

};

//=============================================================

// Rule with id
class Rule_id: public Rule {
public:
    Rule_id(Grammar* gr, const string& s, int i): Rule(s) {
        id =i;
        name = s+"["+to_string(id)+"]";
        grammar = gr;
        grammar->id_names[id] = name;
    }

    virtual Rule_result do_match(bool do_ast) override {
        auto res =  rules.front()->match(do_ast);
        // If Ok continue with ast
        if(res.ok){
            if(do_ast){
                if(!res.ast){
                    // If no ast returned create ast node and put matched string inside
                    res.ast.reset(new Ast_node(id));
                    res.ast->value = string(res.first,res.last);
                } else {
                    // Set id for returned ast
                    res.ast->id = id;
                    // Flatten anonimous nodes
                    if(!res.ast->nodes.empty()){
                        auto it=res.ast->nodes.begin();
                        do {
                            cout << "::" << (*it)->id << endl;
                            if((*it)->id<0){
                                res.ast->nodes.insert(it,(*it)->nodes.begin(),(*it)->nodes.end());
                                res.ast->nodes.erase(it);
                                it=res.ast->nodes.begin(); // restant since list changed
                            }
                            it++;
                        } while(it!=res.ast->nodes.end());
                    }
                }
            }
        }
        return res;
    }
};

// Rule which supresses AST generation recursively
class Rule_no_ast: public Rule {
public:
    Rule_no_ast(const string& s): Rule(s) {}

    virtual Rule_result do_match(bool do_ast) override {
        return rules.front()->match(false);
    }
};


// Matches literaly
class Rule_literal: public Rule {
public:
    Rule_literal(const string& s): Rule("lit:"+s), lit(s) {
        kind=Kind::literal;
    }

    virtual Rule_result do_match(bool do_ast) override {
        string::const_iterator inp_it = grammar->cur;
        string::const_iterator lit_it = lit.begin();
        Rule_result res;

        while(*inp_it==*lit_it){
            inp_it++;
            lit_it++;
            if(lit_it==lit.end()){
                res.ok = true; //success
                res.last = inp_it;
                return res;
            } else if(inp_it==grammar->end_input) {
                res.ok = false; //failed
                res.last = inp_it;
                return res;
            }
        }
        // Failed if here
        res.ok = false;
        res.last = inp_it;
        return res;
    }
private:
    string lit;
};

// Matches any charater
class Rule_any_char: public Rule {
public:
    Rule_any_char(const string& s): Rule("any_ch:"+s), chars(s) {}

    virtual Rule_result do_match(bool do_ast) override {
        Rule_result res;
        res.last = grammar->cur+1;
        for(char c: chars){
            if(c==*grammar->cur){
                res.ok=true;
                return res;
            }
        }
        // Failed if here
        res.ok = false;
        return res;
    }
private:
    string chars;
};



// Matches sequence of subrules
class Rule_sequence: public Rule {
public:
    Rule_sequence(): Rule() {kind=Kind::sequence; name="seq";}

    virtual Rule_result do_match(bool do_ast) override {
        Rule_result res;
        res.ast.reset(new Ast_node(-1)); // Create temporary ast

        // Try to match all subrules
        for(auto& r: rules){
            auto out = r->match(do_ast);
            res.last = out.last;
            if(!out.ok){ // failure
                res.ok = false;
                res.ast.reset();
                return res;
            } else { // success
                if(out.ast) res.ast->nodes.push_back(out.ast);
            }
        }
        // Success if here
        if(res.ast->nodes.empty()){
            res.ast.reset();
        } else if(res.ast->nodes.size()==1) res.ast = res.ast->nodes.front();
        res.ok = true;
        return res;
    }
};


// Matches alternatives
class Rule_alternative: public Rule {
public:
    Rule_alternative(): Rule() {kind=Kind::alternative; name="alt";}

    virtual Rule_result do_match(bool do_ast) override {
        // Try to match
        Rule_result out;
        for(auto& r: rules){
            auto out = r->match(do_ast);
            if(out.ok) return out; // If matched return immediately
        }
        // Failure if here
        return out;
    }
};


// Matches on or more times
class Rule_one_or_more: public Rule {
public:
    Rule_one_or_more(): Rule() {name="one_or_more";}

    virtual Rule_result do_match(bool do_ast) override {
        Rule_result fin; //final result
        Rule_result res = rules.front()->match(do_ast);
        if(!res.ok) return res;

        fin.last = res.last;
        fin.ok = res.ok;
        if(res.ast){
            fin.ast.reset(new Ast_node(-1));
            fin.ast->nodes.push_back(res.ast);
        }

        do {
            res = rules.front()->match(do_ast);
            if(res.ok){
                fin.last = res.last;
                if(res.ast) fin.ast->nodes.push_back(res.ast);
            }
        } while(res.ok);
        return fin;
    }
};

class Rule_zero_or_more: public Rule {
public:
    Rule_zero_or_more(): Rule() {name="zero_or_more";}

    virtual Rule_result do_match(bool do_ast) override {
        Rule_result res = rules.front()->match(do_ast);
        if(!res.ok){ // If not matched it's fine
            res.ok = true;
            return res;
        }
        Rule_result out;
        do {
            out = rules.front()->match(do_ast);
            if(out.ok) res = out;
        } while(out.ok);
        return res;
    }
};

class Rule_zero_or_one: public Rule {
public:
    Rule_zero_or_one(): Rule() {name="zero_or_one";}

    virtual Rule_result do_match(bool do_ast) override {
        auto saved = grammar->cur;
        Rule_result res = rules.front()->match(do_ast);
        // If not matched backtrack here sine ok will be returned any way
        if(!res.ok) res.last = saved;
        res.ok = true; // Ok any way
        return res;
    }
};


/*
class Rule_ast: public Rule {
public:
    Rule_ast(): Rule() {name="ast";}

    virtual Rule_result do_match() override {
        // Pass flag to all tree under this rule
        // If any rule in the tree has
        Rule_result res = rules.front()->match(true);
        // If ok generate ast node
        if(res.ok){
            res.ast
        }
        return res;
    }
};
*/


void Rule::set_grammar(Grammar *gr){
    //grammar = gr;
    //for(auto& r: rules) r->set_grammar(gr);
}

string Rule::to_str(string prefix){
    string res = prefix + name;
    if(enable_ast) res+="!";
    //if(id<0){
    if(!rules.empty()) res+=" {\n";
    for(auto& r: rules) res += r->to_str(prefix+" ");
    if(!rules.empty())
        res+= prefix+"}\n";
    else
        res+="\n";
    //} else res+="\n";
    return res;
}

Rule_result Rule::match(bool gen_ast){
    // Save current parser position
    auto old = grammar->cur;
    // Call matcher

    cout << "try " << name << " at " << *grammar->cur << endl;

    auto res = do_match(gen_ast);
    res.first = old;
    // If failed retract
    if(!res.ok){
        grammar->cur = old;
    } else {
        // success, advance the position
        grammar->cur = res.last;

        cout << "matched "<< name << "(" << string(res.first,res.last) << ")"
             << " ast=" << res.ast << endl;

        // If success and asked to skip spaces consume as many as we can
        if(res.ok && grammar->skip_whitespaces){
            while(isspace(*grammar->cur)) grammar->cur++;
        }
    }
    return res;
}


/*
Rule_holder operator+(const Rule_holder& cur){
    // Enable AST generation for this rule if it doesn't have an id
    if(cur->id<0){
        cur->enable_ast = true;
    }
    return cur;
}
*/


Rule_holder Grammar::lit(const string &s){ return add_rule(new Rule_literal(s)); }

Rule_holder Grammar::ch(const string &s){ return add_rule(new Rule_any_char(s)); }


//========================================================

class MyGrammar: public Grammar {
public:

    RULE(INTEGER) <= ch("+-")['?'] << ch("0123456789")['+'];
    RULE(FLOAT) <= INTEGER << (lit(".") << ch("0123456789")['+'])['?'];
    RULE(SCI) <= !(FLOAT << (ch("eE") << INTEGER)['?']);
    RULE(OP_MUL) <= lit("*");
    RULE(OP_DIV) <= lit("/");
    RULE(OP_PLUS) <= lit("+");
    RULE(OP_MINUS) <= lit("-");
    RULE(FACTOR);
    RULE(TERM);
    RULE(EXPR);

    MyGrammar(): Grammar() {
        FACTOR <= SCI | (lit("(") << EXPR << lit(")"));
        TERM <= FACTOR << ((OP_MUL|OP_DIV) << FACTOR)['*'];
        EXPR <= TERM << ((OP_PLUS|OP_MINUS) << TERM)['*'];

        top_rule = EXPR.get();
    }
};


int main(int argc, char* argv[]){
    MyGrammar g;
    g.set_input("2+2");

    //cout << g.to_string() << endl;
    cout << "Number of rules in grammar:" << g.rules.size() << endl;

    auto res = g.parse();

    cout << res.ok << endl;
    cout << "AST:" << endl;
    cout << res.ast->to_str(&g) << endl;

        //for(list<unique_ptr<Rule>>::iterator it=g.rules.begin(); it!=g.rules.end();it++){
    //    cout << (*it)->to_str() << endl;
    //}
}

Rule_holder Rule_holder::operator<<(const Rule_holder &rhs){
    if(rule->kind==Rule::Kind::sequence){
        rule->rules.push_back(rhs);
        return *this;
    } else {
        // create new parent
        Rule_holder parent = rule->grammar->add_rule(new Rule_sequence());
        parent->rules.push_back(*this);
        parent->rules.push_back(rhs);
        return parent;
    }
}

Rule_holder Rule_holder::operator|(const Rule_holder &rhs){
    if(rule->kind==Rule::Kind::alternative){
        rule->rules.push_back(rhs);
        return *this;
    } else {
        // rule is not alternative
        Rule_holder parent = rule->grammar->add_rule(new Rule_alternative());
        parent->rules.push_back(*this);
        parent->rules.push_back(rhs);
        return parent;
    }
}

Rule_holder Rule_holder::operator[](char c){
    // Match modifier
    switch(c){
    case '+': {
        Rule_holder r = rule->grammar->add_rule(new Rule_one_or_more());
        r->rules.push_back(*this);
        return r;
    }
    case '*': {
        Rule_holder r = rule->grammar->add_rule(new Rule_zero_or_more());
        r->rules.push_back(*this);
        return r;
    }
    case '?': {
        Rule_holder r = rule->grammar->add_rule(new Rule_zero_or_one());
        r->rules.push_back(*this);
        return r;
    }
    }
}

Rule_holder Rule_holder::operator!(){
    auto r = rule->grammar->add_rule(new Rule_no_ast("!"+rule->name));
    r->rules.push_back(*this);
    return r;
}

string Ast_node::to_str(Grammar *gr, string prefix){
    string res = prefix + gr->id_names[id] + ": " + get<string>(value);
    if(!nodes.empty()) res+=" {\n";
    for(auto& n: nodes) res += n->to_str(gr, prefix+" ");
    if(!nodes.empty())
        res+= prefix+"}\n";
    else
        res+="\n";
    //} else res+="\n";
    return res;
}
