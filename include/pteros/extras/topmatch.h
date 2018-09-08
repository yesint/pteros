#pragma once

#include "pteros/core/selection.h"
#include "pteros/core/logging.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>
#include <set>
#include <list>
#include <deque>

namespace pteros {


struct Node {
    Node(){}
    Node(int par): parent(par) {}
    int parent;
    std::vector<int> children;
};

struct Tree {
    Tree(int ind): root_ind(ind) {
        tree[root_ind] = Node(-1);
    }

    void add(int parent, int ind){
        if(tree.count(parent)==0)
            throw Pteros_error("No such parent: {}",parent+1);
        if(tree.count(ind)>0) throw Pteros_error("Index exists: {}",ind+1);
        tree[ind] = Node(parent);
        tree[parent].children.push_back(ind);
    }

    Node& operator[](int i){
        if(tree.count(i)==0) throw Pteros_error("No such index in tree: {}",i+1);
        return tree[i];
    }

    bool has(int i){ return tree.count(i)>0; }

    void print(int cur=-1, std::string pad="");

    std::map<int,Node> tree;
    int root_ind;    
};


//-------------------------------------------

class Topmatch {
public:
    Topmatch(){}
    Topmatch(const Selection& sel);
    void set_source(const Selection& sel);

private:
    void build_tree(int tree_ind, int cur);

    Selection* p_sel;
    std::vector<std::vector<int>> con;

    std::vector<Tree> trees;
    std::deque<Eigen::Vector2i> todo;
};

void use_babel(const Selection& sel);

}
