/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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
            mask(sel.Index(i)) = 1;
        }
        b = sel.Index(0);
        e = sel.Index(sel.size()-1);
        coord = &sel.get_system()->Frame_data(sel.get_frame()).coord;
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
        for(int i=0; i<10; ++i)
            system().rotate(1,0.1);

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
    for(int i=0;i<sel.size();++i) (*data)[i] = sel.XYZ(i);
    return [data](int i) -> Vector3f& {
        return (*data)[i];
    };
}



int main(int argc, char** argv)
{    

    vector<Lipid_descr> species = {
        {"YOPE","resname YOPE", "name P N", "name C314 C315 C316 C216 C217 C218", "name C22 C21 C23 C31 C32 C33"},
        {"DYPE","resname DYPE", "name P N", "name C314 C315 C316 C214 C215 C216", "name C22 C21 C23 C31 C32 C33"},
        {"DOPE","resname DOPE", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33"},
        {"DOPC","resname DOPC", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33"},
        {"POPE","resname POPE", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33"},
        {"DYPC","resname DYPC", "name P N", "name C314 C315 C316 C214 C215 C216", "name C22 C21 C23 C31 C32 C33"},
        {"YOPC","resname YOPC", "name P N", "name C314 C315 C316 C216 C217 C218", "name C22 C21 C23 C31 C32 C33"},
        {"POPC","resname POPC", "name P N", "name C316 C317 C318 C216 C217 C218", "name C22 C21 C23 C31 C32 C33"},
    };

    System sys("/home/semen/work/current/Projects/Masato/symmetric/topol.tpr");
    sys().write("a.pdb");

    vector<Selection> res;
    sys().split_by_molecule(res);
    cout << res.size() << endl;


    /*
    System sys("/home/semen/work/current/Projects/Masato/symmetric/after_em.gro");

    Membrane membr(&sys,species);
    membr.compute_properties(2.0, Vector3f(0,0,1));
    membr.write_vmd_arrows("normals.tcl");
    membr.write_smoothed("smoothed.pdb");    

    // lip properties

    for(int i=0;i<membr.lipids.size();++i){
        cout << i << " " << membr.lipids[i].normal(2) << " " << membr.lipids[i].area << endl;
    }
    */

    return 1;

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
        for(int i=0;i<sel.size();++i) data[i] = sel.XYZ(i);


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
        for(int i=0;i<sel.size();++i) data.col(i) = sel.XYZ(i);

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
                sel.XYZ(i)=sel.XYZ(i-1)+sel.XYZ(j);
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
                sel[i].XYZ()=sel[i-1].XYZ()+sel[j].XYZ();
            }
        }

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout << "Execution time: " << duration << endl;
    }
    //-----------------------

    return 1;
*/
    try{



        Options opt;
        parse_command_line(argc,argv,opt);
        Trajectory_reader reader(opt);

        reader.add_task( new Test_serial(opt) );
        reader.add_task( new Test_serial(opt) );

        //reader.add_task( new Test_serial(opt) );
        //reader.add_task( new Test_serial(opt) );

          //reader.add_task( new Test_serial() );
        //reader.add_task( new Test_serial() );

        reader.run();


        //cout << (boost::get<Parse_tree_ptr>(p->children.front())) << endl;

        //std::shared_ptr<Parser> p(new Parser);
        //p->parse();

//        System s("/home/semen/work/Projects/asymmetric_bilayer/for-diamonds/hexagonal/2x2.gro");
  //      Selection sel(s,"(y+4)<(x+2) or (1-z)>x");
/*
        System s("/home/semen/work/Projects/pteros/paper2/supplement/benchmark/large.gro");


        int Ntry = 1000;


        Vector3f r1(Vector3f::Zero()), r2(Vector3f::Zero());

        ofstream f("size-dep.dat");

        for(int N=2;N<s.num_atoms();N+=50000){
            // Make sparse vector
            vector<int> ind;
            ind.reserve(N);
            for(int i=0;i<N;++i){
                ind.push_back(i*(s.num_atoms()/(N-1)));
            }

            // Normal selection
            r1.fill(0);
            Selection sel(s,ind);

            auto t_start = std::chrono::high_resolution_clock::now();
            for(int n=0;n<Ntry;++n)
                r1 += sel.center();
            auto t_end = std::chrono::high_resolution_clock::now();
            auto T1 = 1e6*std::chrono::duration<double>(t_end-t_start).count()/float(Ntry);

            // Mask access
            r2.fill(0);
            SelTest_mask sel1(sel);

            t_start = std::chrono::high_resolution_clock::now();
            for(int n=0;n<Ntry;++n)
                r2 += sel1.center();
            t_end = std::chrono::high_resolution_clock::now();
            auto T2 = 1e6*std::chrono::duration<double>(t_end-t_start).count()/float(Ntry);

            cout << N << " " << T1 << " " << T2 << " | " << r1.transpose() << " :: " << r2.transpose() << endl;
            f << N << " " << T1 << " " << T2 << endl;
        }

        f.close();
*/

/*
        Selection w;
        Selection sel(s,"all");
        int N = 10;
        vector<Vector2i> bon;
        float d = 0.6;

        auto t_start = std::chrono::high_resolution_clock::now();
        for(int i=0;i<N;++i)
            w.modify(s,"within 0.6 nopbc of name PO4");
            //Grid_searcher(d,sel,bon);
        auto t_end = std::chrono::high_resolution_clock::now();

        cout << " elapsed: "
             << 1e6*std::chrono::duration<double>(t_end-t_start).count()/float(N) << endl;
*/

    } catch(const Pteros_error& e){
        LOG()->error(e.what());
    }

}

