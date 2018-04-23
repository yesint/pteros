#pragma once

#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>
#include <set>
#include <list>

namespace pteros {

struct Mol_node {
    Mol_node(){}
    Mol_node(int _par, int _ind, int _el): parent(_par), ind(_ind), element(_el) { }
    std::list<Mol_node>* add_variant(){
        children.emplace_back();
        return &children.back();
    }

    Mol_node* add_child(int n, int el, std::list<Mol_node>* var){
        var->emplace_back(ind,n,el);
        return &var->back();
    }

    void print(int tab=0);
    void get_ind_vector(std::vector<int>& v);

    int ind;
    int element;
    int parent;
    std::list<std::list<Mol_node>> children;
};

//-------------------------------------------

class Topmatch {
public:
    Topmatch(){}
    Topmatch(const Selection& sel);
    void set_source(const Selection& sel);
    // Match another molecule against this one
    bool match(const Selection& sel);
    // Match molecule against itself to determine symmetry
    int match_self();

    std::vector<int> get_mapping();

private:
    void build_tree(Mol_node& node, std::shared_ptr<std::set<int>>& used);
    bool build_match(Mol_node& node, const Mol_node& ref);

    std::vector<std::vector<int>> con, m_con;
    Selection* p_sel;    
    Mol_node root, m_root;
};

}
