
#include "pteros/pteros.h"
#include "pteros/python/compiled_plugin.h"
#include <fstream>

using namespace std;
using namespace pteros;
using namespace Eigen;


TASK_SERIAL(charge_distrib)
public:

    string help() override {
        return
R"(Purpose:
    Compute 1D or 2D charge distribution across the box.
    Only rectangular boxes are supported.
Options:
    -sel <sel_str1 sel_str2 ...>
        Selections to compute charge distributions for.
        Dynamic selections are allowed.
    -Ngrid <n>
        Number of bins across each dimension.
    -oneD <true|false>
        Compute 1D charge distribution.
    -axis <n>
        If oneD is true computes charge distribution across
        this dimension.
        Otherwise indicates a normal dimension to the plane
        where distribution have to be computed.
Output:
    If oneD is true writes a file q_distrib_1D_<n>.dat
    for each selection (n starts from zero).
    This file has columns (x,q), where x is coordinate
    and q is an average charge in atomic units.

    If oneD is false writes a file q_distrib_<n>.dat
    for each selection (n starts from zero).
    This file contans a (Ngrid x Ngrid) matrix containing
    average charge in atomic units.
    Also writes a file q_distrib_extents_<n>.dat
    for each selection, which contains average box sizes
    in the plane where average charge is computed.
)";
    }

protected:

    void pre_process() override {
        if(system.box(0).is_triclinic()) throw PterosError("Only rectangular boxes are allowed!");    
    
        auto sel_texts = options("sel").as_strings();
        for(const auto& s: sel_texts) sel.emplace_back(system,s);
        
        _1D = options("oneD","false").as_bool();

        Np = options("Ngrid","100").as_int();
        if(_1D){
            for(int i=0;i<sel.size();++i) q_distrib_1D.emplace_back(0,1,Np);
        } else {
            for(int i=0;i<sel.size();++i) q_distrib.emplace_back(0,1,Np,0,1,Np);
        }
        
        axis = options("axis","2").as_int();
        
        if(!_1D){
            // Set working dimensions
            switch(axis){
                case 0:
                    ind1 = 1; ind2 = 2;
                    break;
                case 1:
                    ind1 = 0; ind2 = 2;
                    break;
                case 2:
                    ind1 = 0; ind2 = 1;
                    break;
            }
        } //if    
        
        ave_box_sz.fill(0.0);    
    }


    void process_frame(const FrameInfo &info) override {
        for(auto& s: sel){
            s.apply();
            s.wrap();
        }
                       
        auto ext = system.box().extents();
        
        if(_1D){
            float slab_v = system.box().volume()/float(Np);
            
            for(int j=0;j<sel.size();++j){
                // Loop over selection
                for(int i=0;i<sel[j].size();++i){            
                    q_distrib_1D[j].add(sel[j].xyz(i)(axis)/ext(axis), sel[j].charge(i)/slab_v);
                }
            }
            ave_box_sz(0) += ext(axis);
        } else { // 2D
            float cell_v = system.box().volume()/float(Np*Np);           
            
            for(int j=0;j<sel.size();++j){
                // Loop over selection
                for(int i=0;i<sel[j].size();++i){
                    const auto& coor = sel[j].xyz(i);
                    q_distrib[j].add(coor(ind1)/ext(ind1), coor(ind2)/ext(ind2), sel[j].charge(i)/cell_v);
                }
            }
            ave_box_sz(0) += ext(ind1);
            ave_box_sz(1) += ext(ind2);
        } // if 1D
    }


    void post_process(const FrameInfo &info) override {
        //const float EPS0=8.85419E-12;
        //const float ELC=1.60219E-19;
        //const float Fcoeff = ELC * -1.0E9 / EPS0;
        if(_1D){
            ave_box_sz /= float(info.valid_frame);
            for(int j=0;j<sel.size();++j){
                q_distrib_1D[j].normalize(info.valid_frame);
                q_distrib_1D[j].positions() *= ave_box_sz(0);
                q_distrib_1D[j].save_to_file(fmt::format("q_distrib_1D_{}.dat",j));
            }            
        } else { // 2D
            for(int j=0;j<sel.size();++j){
                q_distrib[j].normalize(info.valid_frame);
                q_distrib[j].save_to_file(fmt::format("q_distrib_{}.dat",j));
                // File with average box extents
                ofstream out(fmt::format("q_distrib_extents_{}.dat",j));
                fmt::print(out,"{}\n{}\n", 
                    ave_box_sz(0)/float(info.valid_frame),
                    ave_box_sz(1)/float(info.valid_frame)
                    );
                out.close();
            }
        }
    }


private:
    vector<Histogram2D> q_distrib;
    vector<Histogram> q_distrib_1D;
    vector<Selection> sel;
    int Np,axis;
    bool _1D;
    Vector2f ave_box_sz;
    int ind1,ind2;
};

CREATE_COMPILED_PLUGIN(charge_distrib)

