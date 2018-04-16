#include "pteros/extras/topmatch.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <set>
#include <queue>

using namespace std;
using namespace pteros;
using namespace Eigen;

struct Mol_node {
    Mol_node(int _par, int _ind, int _el): parent(_par), ind(_ind), element(_el) {}

    shared_ptr<Mol_node> add(int i, int e){
        auto el = make_shared<Mol_node>(ind,i,e);
        children.push_back(el);
        return el;
    }

    void print(int tab=0){
        cout << ind+1 << " ";
        //cout << element << " ";
        if(!children.empty()){
            cout << "(";
            for(auto c: children){
                c->print(tab+2);
            }
            cout << ")";
        }
        if(tab==0) cout << endl;
    }

    void get_ind_vector(vector<int>& v){
        v.push_back(ind);
        for(auto c: children) c->get_ind_vector(v);
    }

    int ind;
    int element;
    int parent;
    vector<shared_ptr<Mol_node>> children;
};

typedef shared_ptr<Mol_node> Mol_node_ptr;


//---------------------------

class Tree_builder {
public:
    Tree_builder(const Selection& sel) {
        p_sel = const_cast<Selection*>(&sel);

        vector<Vector2i> bonds;
        search_contacts(0.18,sel,bonds,false,false); // local indexes
        con.resize(sel.size());
        for(auto& b: bonds){
            con[b[0]].push_back(b[1]);
            con[b[1]].push_back(b[0]);
        }

        // Sort each entry by number of bonds
        for(auto& el: con){
            sort(el.begin(),el.end(), [this](int a,int b){return con[a].size()<con[b].size();});
        }

        // Create root
        root = make_shared<Mol_node>(-1,0,sel.element_number(0));
        build_tree(root);
    }

    Mol_node_ptr get_tree(){ return root; }

    // Match another molecule against this one
    bool match(const Selection& sel){
        p_sel = const_cast<Selection*>(&sel);

        vector<Vector2i> bonds;
        search_contacts(0.18,sel,bonds,false,false); // local indexes

        m_con.clear();
        m_con.resize(sel.size());
        for(auto& b: bonds){
            m_con[b[0]].push_back(b[1]);
            m_con[b[1]].push_back(b[0]);
        }

        // Sort each entry by number of bonds
        for(auto& el: m_con){
            sort(el.begin(),el.end(), [this](int a,int b){return m_con[a].size()<m_con[b].size();});
        }

        for(int i=0; i<sel.size(); ++i){
            used.clear();
            m_root = make_shared<Mol_node>(-1,i,sel.element_number(i));                        
            if(build_match(m_root,root)) return true;
        }
        // If we are here no match is found
        return false;
    }

    vector<int> get_mapping(){
        vector<int> ind1, ind2, res;

        root->get_ind_vector(ind1);
        m_root->get_ind_vector(ind2);

        // Get sort permutation for ind1
        vector<std::size_t> per(ind1.size());
        iota(per.begin(),per.end(), 0);
        sort(per.begin(),per.end(), [&ind1](int i,int j){return ind1[i]<ind1[j];});
        // Apply permutation to ind2 in place
        res.resize(ind1.size());
        for(int i=0;i<ind1.size();++i){
            res[i] = ind2[per[i]];
        }

        return(res);
    }


private:
    void build_tree(Mol_node_ptr& node){
        used.insert(node->ind);

        for(int b: con[node->ind]){
            if(node->parent==b) continue;
            if(used.count(b)) continue;
            auto el = node->add(b,p_sel->element_number(b));
            build_tree(el);
        }
    }



    bool build_match(Mol_node_ptr& node, Mol_node_ptr& ref){
        if(node->element != ref->element) return false;
        if(m_con[node->ind].size() != con[ref->ind].size()) return false;

        used.insert(node->ind);

        for(int b: m_con[node->ind]){
            if(node->parent==b) continue;
            if(used.count(b)) continue;

            auto el = node->add(b,p_sel->element_number(b));

            bool ok;
            for(auto& r: ref->children){
                ok = build_match(el,r);
                if(ok) break;
            }
            if(!ok){
                used.erase(node->ind);
                node->children.clear();
                return false;
            }
        }
        return true; // All children matched
    }

    vector<vector<int>> con, m_con;
    Selection* p_sel;
    set<int> used;
    Mol_node_ptr root, m_root;
};

//---------------------------


void pteros::Topmatch::compute_mapping(const pteros::Selection &src, const pteros::Selection &target)
{    
    // sanity check
    if(src.size()!=target.size()) throw Pteros_error("Can't match selections of different size");

    vector<int> mapping;

    Tree_builder builder(src);
    if(builder.match(target)){
        mapping = builder.get_mapping();

        for(auto el: mapping) cout << el << endl;
    }
}
