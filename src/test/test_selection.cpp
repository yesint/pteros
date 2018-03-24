/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *  
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#include <string>
#include "pteros/analysis/options.h"
#include <Eigen/Core>
#include "pteros/core/pteros_error.h"
#include "pteros/core/selection.h"
#include "pteros/core/distance_search.h"


#include "pteros/analysis/trajectory_reader.h"
#include "pteros/analysis/task_plugin.h"

#include <chrono>
#include <fstream>
#include <thread>

#include "pteros/core/logging.h"
#include "spdlog/fmt/ostr.h"

#include "pteros/extras/membrane.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

//-------------------------------------------------------


struct Sel_expr {
    string name;
    vector<Sel_expr*> children;
    virtual string dump(int level=0){
        string shift;
        for(int i=0;i<level;++i) shift+="  ";

        string s = shift + name + " {\n";
        for(auto& c: children){
            s += c->dump(level+1);
        }
        s += shift+"}\n";
        return s;
    }
    virtual ~Sel_expr(){
        cout << "del: " << name << endl;
        for(int i=0;i<children.size();++i) delete children[i];
    }
};

struct Sel_container {
    Sel_expr* expr;
    Sel_container(Sel_expr* ex): expr(ex){}
    Sel_expr* operator()() const {return expr;}
};


struct Sel_expr_name: public Sel_expr {
    Sel_expr_name(string val): value(val) {
        name = "name";
    }

    virtual string dump(int level=0){
        string shift;
        for(int i=0;i<level;++i) shift+="  ";
        return shift+name+":"+value+"\n";
    }

    string value;
};

Sel_container name(string val){
    return Sel_container(new Sel_expr_name(val));
}

struct Sel_expr_and: public Sel_expr {
    Sel_expr_and(Sel_expr* op1, Sel_expr* op2){
        name = "and";
        children.push_back(op1);
        children.push_back(op2);
    }
};

struct Sel_expr_or: public Sel_expr {
    Sel_expr_or(Sel_expr* op1, Sel_expr* op2){
        name = "or";
        children.push_back(op1);
        children.push_back(op2);
    }
};

struct Sel_expr_not: public Sel_expr {
    Sel_expr_not(Sel_expr* op1){
        name = "not";
        children.push_back(op1);
    }
};

Sel_container operator&(const Sel_container& op1, const Sel_container& op2){
    return Sel_container(new Sel_expr_and(op1(),op2()));
}

Sel_container operator|(const Sel_container& op1, const Sel_container& op2){
    return Sel_container(new Sel_expr_or(op1(),op2()));
}

Sel_container operator!(const Sel_container& op1){
    return Sel_container(new Sel_expr_not(op1()));
}


//-------------------------------------------------------

class SelTest_mask {
public:
    SelTest_mask(const Selection& sel){
        mask.resize(sel.get_system()->num_atoms());
        mask.fill(0);
        n = sel.size();
        for(int i=0;i<sel.size();++i){
            mask(sel.index(i)) = 1;
        }
        b = sel.index(0);
        e = sel.index(sel.size()-1);
        coord = &sel.get_system()->frame(sel.get_frame()).coord;
    }

    ~SelTest_mask(){}

    Vector3f center(){
        Vector3f res(Vector3f::Zero());
        for(int i=b; i<=e; ++i)
            if(mask(i)==1)
                res += (*coord)[i];
        return res/float(n);
    }

private:
    VectorXi mask;
    int b,e,n;
    vector<Vector3f>* coord;
};


//-------------------------------------------------------


TASK_PARALLEL(Test_task)
    void pre_process() override {
        jump_remover.add_atoms(system("index 1-10"));
        log->info("Test_task pre_process #{}", get_id());
        LOG()->info("Generic message");
    }

    void process_frame(const Frame_info& info) override {


        log->info("Test_task process_frame {} #{}", info.valid_frame,get_id());
        //fflush(stdout);

        //std::this_thread::sleep_for(std::chrono::seconds(1));
        //for(int i=0; i<10; ++i)
        //    system().rotate(1,0.1);

        res.push_back(info.valid_frame);
        //throw Pteros_error("UPS");
    }

    void post_process(const Frame_info& info) override {
        log->info("Test_task post_process of instance {}", info.last_frame);
    }

    void collect_data(const vector<Task_ptr>& tasks) override {
        for(auto& t: tasks){
            auto h = dynamic_cast<Test_task*>(t.get());
            Eigen::Map<VectorXi> m(h->res.data(),h->res.size());
            LOG()->info("{} :: {}", h->get_id(), m.transpose());
        }
    }

public:
    vector<int> res;
};


void accum(const Frame_info& info, const std::vector<Task_ptr>& tasks){    
    for(auto& t: tasks){
        auto h = dynamic_cast<Test_task*>(t.get());
        Eigen::Map<VectorXi> m(h->res.data(),h->res.size());
        LOG()->info("{} :: {}", h->get_id(), m.transpose());
    }
}


TASK_SERIAL(Test_serial)
    void pre_process() override {
        cout << "Test_serial pre_process" << endl;
    }

    void process_frame(const Frame_info& info) override {

        cout << "Test_serial process_frame " << std::this_thread::get_id() << " "<< info.valid_frame << endl;
        //std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    void post_process(const Frame_info& info) override {
        cout << "Test_serial post_process" << info.last_frame << endl;
    }
};

//-----------------------------

std::function<Vector3f&(int)>
make_accessor(const Selection& sel){
    shared_ptr<vector<Vector3f>> data = make_shared<vector<Vector3f>>(sel.size());
    for(int i=0;i<sel.size();++i) (*data)[i] = sel.xyz(i);
    return [data](int i) -> Vector3f& {
        return (*data)[i];
    };
}
//-----------------------------

int main(int argc, char** argv)
{    

    /*
    string t1("name C21 C22 C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216");
    string t2("name C21 C22 C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216");

    vector<Lipid_descr> species = {
        {"YOPE","resname YOPE", "name P N", "name C314 C315 C316 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"DYPE","resname DYPE", "name P N", "name C314 C315 C316 C214 C215 C216", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"DOPE","resname DOPE", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"DOPC","resname DOPC", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"POPE","resname POPE", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"DYPC","resname DYPC", "name P N", "name C314 C315 C316 C214 C215 C216", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"YOPC","resname YOPC", "name P N", "name C314 C315 C316 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
        {"POPC","resname POPC", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33", {t1,t2}},
    };

    System sys("/home/semen/work/current/Projects/Masato/symmetric/after_em.gro");

    Membrane membr(&sys,species);
    membr.compute_properties(2.0, true, Vector3f(0,0,0), Vector3i(0,0,1));
    membr.write_vmd_arrows("normals.tcl");
    membr.write_smoothed("smoothed.pdb");    

    // lip properties

    for(int i=0;i<membr.lipids.size();++i){
        //cout << Map<VectorXf>(membr.lipids[i].order[0].data(),16).transpose() << endl;
        cout << membr.lipids[i].area << endl;
    }

    return 1;
    */

    //System s("/home/semen/work/current/Projects/plasma_vesicle/after_em.gro");
/*
    // System with 10^6 atoms
    System s("/home/semen/work/stored/Projects/pteros/paper2/supplement/benchmark/large.gro");
    vector<int> ind;
    for(int i=0;i<s.num_atoms();i+=1) ind.push_back(i);
    Selection sel(s,ind);

    std::clock_t start;
    double duration;

    //-----------------------

    {
        start = std::clock();

        vector<Vector3f> data(sel.size());
        for(int i=0;i<sel.size();++i) data[i] = sel.xyz(i);


        for(int i=1;i<sel.size();++i){
            for(int j=1;j<1000;++j){
                data[i]=data[i-1]+data[j];
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;
    }

    //-----------------------

    {
        start = std::clock();

        MatrixXf data(3,sel.size());
        for(int i=0;i<sel.size();++i) data.col(i) = sel.xyz(i);

        for(int i=1;i<sel.size();++i){
            for(int j=1;j<1000;++j){
                data.col(i)=data.col(i-1)+data.col(j);
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;
    }

    //-----------------------
    {

        start = std::clock();

        for(int i=1;i<sel.size();++i){
            for(int j=1;j<1000;++j){
                sel.xyz(i)=sel.xyz(i-1)+sel.xyz(j);
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;
    }

    //-----------------------


    {

        start = std::clock();

        for(int i=1;i<sel.size();++i){
            for(int j=1;j<1000;++j){
                sel[i].xyz()=sel[i-1].xyz()+sel[j].xyz();
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;
    }
    //-----------------------

    return 1;
*/
    try{


        System s("/home/semen/work/current/Projects/Ache/cc_prot/topol.tpr");
        for(auto& a: s()) cout << a.vdw() << endl;


    } catch(const Pteros_error& e){
        LOG()->error(e.what());
    }

}


