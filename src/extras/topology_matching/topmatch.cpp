#include "pteros/extras/topmatch.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <set>
#include <queue>

/*
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace pteros;
using namespace Eigen;
using namespace OpenBabel;

OBMol sel_to_obmol(const Selection& sel){
    OBMol mol;

    mol.BeginModify();

    for(int i=0;i<sel.size();++i){
        auto& at = sel.atom(i);
        // Create new atom in this mol
        auto oba = mol.NewAtom();
        oba->SetAtomicNum(at.element_number);
        //oba->SetPartialCharge(at.charge);
        oba->SetVector(sel.x(i),sel.y(i),sel.z(i));
    }

    //mol.ConnectTheDots();
    //mol.PerceiveBondOrders();

    mol.EndModify();    

    return mol;
}


void pteros::Topmatch::compute_mapping(const pteros::Selection &src, const pteros::Selection &target)
{
    // sanity check
    if(src.size()!=target.size()) throw Pteros_error("Can't match selections of different size");
    //OBConversion conv;
    //OBMol src_mol, target_mol;
    //conv.ReadFile(&src_mol,"/home/semen/work/current/Projects/Ache/mhc-opt.mol2");
    //conv.ReadFile(&target_mol,"/home/semen/work/current/Projects/Ache/mhc-opt.mol2");

    auto src_mol = sel_to_obmol(src);
    auto target_mol = sel_to_obmol(target);
    FOR_ATOMS_OF_MOL(a,src_mol)
    {
        cout << a->GetX() << endl;
    }


    OBQuery *query = CompileMoleculeQuery(&src_mol);
    cout << "(1)" << endl;
    OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
    cout << "(2)" << endl;
    OBIsomorphismMapper::Mapping mapping;
    mapper->MapFirst(&target_mol,mapping);
    cout << "(3)" << endl;

    for(auto& el: mapping){
        cout << el.first << " " << el.second << endl;
    }

    delete mapper;
    delete query;
}
*/

using namespace std;
using namespace pteros;
using namespace Eigen;

struct Mol_node {
    Mol_node(int _ind, int _el): ind(_ind), element(_el) {}

    shared_ptr<Mol_node> add(int i, int e){
        auto it = make_shared<Mol_node>(i,e);
        children.push_back(it);
        return it;
    }

    void print(int tab=0){
        if(!children.empty()){
            cout << element << "(";
            for(auto c: children){
                c->print(tab+2);
            }
            cout << ")";
        } else {
            cout << element << " ";
        }
        if(tab==0) cout << endl;
    }

    int ind;
    int element;
    vector<shared_ptr<Mol_node>> children;
};

typedef shared_ptr<Mol_node> Mol_node_ptr;


shared_ptr<Mol_node> build_tree(const vector<Vector2i>& bonds, const Selection& sel){
    // Convert bonds toconvinient struct
    vector<vector<int>> con(sel.size());
    for(auto& b: bonds){
        con[b[0]].push_back(b[1]);
        con[b[1]].push_back(b[0]);
    }

    // Create root
    auto root = make_shared<Mol_node>(0,sel.element_number(0));

    set<int> used;

    queue<shared_ptr<Mol_node>> todo;
    todo.push(root);

    while(!todo.empty()){
        auto cur = todo.front();
        todo.pop();
        used.insert(cur->ind);
        // Add bonds to current node
        for(int b: con[cur->ind]){
            if(used.count(b)==0){
                auto el = cur->add(b,sel.element_number(b));
                todo.push(el);
                //used.insert(b);
            }
        }
    }

    return root;
}




bool match_nodes(shared_ptr<Mol_node>& cur, const shared_ptr<Mol_node>& ref){
    if(cur->element != ref->element){
        return false;
    } else if(cur->children.size() != ref->children.size()){
        return false;
    } else if(cur->children.size()==0){
        // For equal lefs
        return true;
    } else {
        // For trees use recursion
        vector<shared_ptr<Mol_node>> matched(cur->children.size());

        vector<int> used(cur->children.size());
        for(int i=0; i<ref->children.size(); ++i) used[i]=0;

        int num = 0;
        // Try to match each child of cur with all unused children of ref
        for(auto& c: cur->children){
            for(int i=0; i<ref->children.size(); ++i){
                if(!used[i] && match_nodes(c,ref->children[i])){
                    matched[i] = c;
                    used[i] = 1;
                    ++num;
                }
            }
        }

        if(num == ref->children.size()){
            // All children matched!
            // Set children in matched order
            cur->children = matched;
            return true;
        } else {
            return false;
        }
    }
}

//---------------------------

class Tree_builder {
public:
    Tree_builder( Selection& sel): p_sel(&sel) {
        vector<Vector2i> bonds;
        search_contacts(0.16,sel,bonds,false,false); // local indexes
        con.resize(sel.size());
        for(auto& b: bonds){
            con[b[0]].push_back(b[1]);
            con[b[1]].push_back(b[0]);
        }

        // Create root
        auto root = make_shared<Mol_node>(0,sel.element_number(0));
        used.insert(0);
        for(int b: con[0]){
            add_child(root,b);
        }
    }

    void add_child(Mol_node_ptr node, int ind){
        if(used.count(ind)>0) return;

        auto el = node->add(ind,p_sel->element_number(ind));
        used.insert(ind);

        for(int b: con[ind]){
            add_child(el,b);
        }
    }

private:
    vector<vector<int>> con;
    Selection* p_sel;
    set<int> used;
};

//---------------------------


void pteros::Topmatch::compute_mapping(const pteros::Selection &src, const pteros::Selection &target)
{
    cout << src.size() << endl;
    // sanity check
    if(src.size()!=target.size()) throw Pteros_error("Can't match selections of different size");

    vector<Vector2i> src_bonds, target_bonds;
    search_contacts(0.16,src,src_bonds,false,false); // local indexes
    search_contacts(0.16,target,target_bonds,false,false); // local indexes

    auto src_tree = build_tree(src_bonds,src);
    auto target_tree = build_tree(target_bonds,target);



    src_tree->print();
}
