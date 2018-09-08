#include "pteros/extras/topmatch.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <set>
#include <numeric>
#include <queue>
#include <numeric>

#include "openbabel/query.h"
#include <openbabel/isomorphism.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace pteros;
using namespace Eigen;


/*
void Mol_node::print(int tab){

    string s = ""; //for(int i=0;i<tab;++i) s+=".";

    cout << s << ind+1 << " ";
    //cout << element << " ";
    int i=0;
    for(auto& var: variants){
        if(variants.size()>1) cout << s << "var " << i << "[" << endl << s;
        if(!var.empty()){
            cout << "(";
            for(auto c: var){
                c.print(tab+2);
            }
            cout << s << ")";
        }
        if(tab==0) cout << endl;
        if(variants.size()>1) cout << endl << s << "] var " << i << endl;
        ++i;
    }
}
*/


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

    // Sort each entry by index
    for(auto& el: con){
        //sort(el.begin(),el.end(), [this](int a,int b){return con[a].size()<con[b].size();});
        sort(el.begin(),el.end());
    }

    // Create root
    trees.emplace_back(0);


    todo.push_front({0,0});

    while(!todo.empty()){
        auto item = todo.front();
        todo.pop_front();
        build_tree(item(0),item(1));
    }


    cout << "# trees = " << trees.size() << endl;
    for(auto& t: trees){
        //if(t.tree.size()==p_sel->size()){
        cout << "--------------------" << endl;
        t.print();
        //}
    }
}


void Topmatch::build_tree(int tree_ind, int cur){    

    // Find valid children
    vector<int> valid;
    for(int b: con[cur]){
        if(trees[tree_ind][cur].parent==b) continue;
        if(trees[tree_ind].has(b)) continue;
        valid.push_back(b);
    }

//    cout << tree_ind << "|" << cur+1 << endl;
/*
        for(int b: valid){
            trees[tree_ind].add(cur,b);
        }
        for(int b: valid){
            build_tree(tree_ind,b);
        }
 */


    if(valid.size()==1){
        trees[tree_ind].add(cur,valid[0]);
        todo.push_front({tree_ind,valid[0]});
    } else if(valid.size()>1){

        auto old = trees[tree_ind];

        int nperm=0;
        do {
            int N;
            if(nperm>0){
                trees.push_back(old);
                N = trees.size()-1;
            } else {
                N = tree_ind;
            }

            for(int b: valid){
                trees[N].add(cur,b);
                todo.push_front({N,b});
            }

            // All todo items already assigned to tree_ind have to be copied to
            // all clones
            if(nperm>0)
                for(auto item: todo){
                    if(item(0)==tree_ind) todo.push_front({N,item(1)});
                }

            ++nperm;
        } while( next_permutation(valid.begin(),valid.end()) );
    }
}

void Tree::print(int cur, string pad)
{
    if(cur<0) cur = root_ind;
    cout << pad << cur+1 << " ";
    if(tree[cur].children.size()>0){
        //cout << "{" << endl;
        cout << endl;
        for(int i=0;i<tree[cur].children.size();++i){
            print(tree[cur].children[i],pad+" ");
        }
        //cout << pad << "}" << endl;
    } else {
        cout << endl;
    }
}

void pteros::use_babel(const Selection &sel)
{
    using namespace OpenBabel;

    OBMol mol;
    sel.to_obmol(mol);

    //OBQuery* query;
    //query = CompileMoleculeQuery(&mol);

    //OpenBabel::OBConversion conv;
    //conv.ReadFile(&mol,"/home/semen/work/current/Projects/Ache/1.pdb");


    std::vector<OBIsomorphismMapper::Mapping> aut;
    FindAutomorphisms(&mol,aut);

    map<int,set<int>> sym;
    for(int i=0;i<sel.size();++i) sym[i]={};

    for(int i=0;i<aut.size();++i){
        cout << "m:"<<i<<endl;
        for(int j=0;j<aut[i].size();++j){
            cout << " " <<  aut[i][j].first+1 << ":" << aut[i][j].second+1 << endl;
            sym[aut[i][j].first].insert(aut[i][j].second);            
        }
    }

    for(auto& it: sym){
        for(int a: it.second){
            if(it.first!=a) sym.erase(a);
        }
    }



    cout << aut.size() << endl;

    for(auto it: sym){
        cout << fmt::format("{} ({}): ",it.first+1,it.second.size());
        for(auto a: it.second) cout << a+1 << " ";
        cout << endl;
    }

    //delete query;
}
