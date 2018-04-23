#include "pteros/extras/topmatch.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <set>
#include <queue>

using namespace std;
using namespace pteros;
using namespace Eigen;



void Mol_node::print(int tab){

    cout << ind+1 << " ";
    //cout << element << " ";
    int i=0;
    for(auto& var: children){
        cout << "var" << i << "[" << endl;
        if(!var.empty()){
            cout << "(";
            for(auto c: var){
                c.print(tab+2);
            }
            cout << ")";
        }
        if(tab==0) cout << endl;
        cout << "]" << endl;
    }



}

void Mol_node::get_ind_vector(vector<int> &v){
    /*
    v.push_back(ind);
    for(auto c: children[0]) c.get_ind_vector(v);
    */
}

Topmatch::Topmatch(const Selection &sel) {
    set_source(sel);
}

void Topmatch::set_source(const Selection &sel)
{
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
    root = Mol_node(-1,0,sel.element_number(0));
    auto used = make_shared<std::set<int>>();
    build_tree(root,used);
    root.print();
}

bool Topmatch::match(const Selection &sel){
    /*
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
        m_root = Mol_node(-1,i,sel.element_number(i));
        if(build_match(m_root,root)) return true;
    }

    // If we are here no match is found
    return false;
    */
}


int Topmatch::match_self(){

}


vector<int> Topmatch::get_mapping(){
    vector<int> ind1, ind2, res;

    root.get_ind_vector(ind1);
    m_root.get_ind_vector(ind2);

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

void Topmatch::build_tree(Mol_node &node, std::shared_ptr<std::set<int>> &used){
    used->insert(node.ind);

    // Find valid children
    vector<int> valid;
    for(int b: con[node.ind]){
        if(node.parent==b) continue;
        if(used->count(b)) continue;
        valid.push_back(b);
    }

    if(valid.size()==1){
        auto var = node.add_variant();
        auto el = node.add_child(valid[0],p_sel->element_number(valid[0]),var);
        build_tree(*el,used); // same used is passed, no copy is made
    } else {

        // Generate all permutations of valid
        do {
            // Create new copy of used and pass it
            auto used_copy = make_shared<set<int>>();
            *used_copy = *used;

            auto var = node.add_variant();
            for(int b: valid){
                auto el = node.add_child(b,p_sel->element_number(b),var);
                build_tree(*el,used_copy);
            }
        } while( next_permutation(valid.begin(),valid.end()) );

    }
}

bool Topmatch::build_match(Mol_node &node, const Mol_node &ref){
    /*
    if(node.element != ref.element) return false;
    if(m_con[node.ind].size() != con[ref.ind].size()) return false;

    used.insert(node.ind);

    for(int b: m_con[node.ind]){
        if(node.parent==b) continue;
        if(used.count(b)) continue;

        auto el = node.add_child(b,p_sel->element_number(b));

        bool ok;
        for(auto& r: ref.children){
            ok = build_match(*el,r);
            if(ok) break;
        }
        if(!ok){
            used.erase(node.ind);
            node.children.clear();
            return false;
        }
    }

    return true; // All children matched
    */
}

