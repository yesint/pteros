/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
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


#include "pteros/python/compiled_plugin.h"
#include <fstream>
#include <map>
#include "pteros/extras/voronoi_packing.h"
#include "pteros/core/system.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


TASK_PARALLEL(voronoi_par)
public:

    string help() override {
        return
R"(Purpose:
    Computes voronoi tesselation for given selections.
Output:
    File voronoi_<id>_groups.dat repors per-group averages
    File voronoi_<id>_interfaces.dat is a matrix of group-group interface areas.
         Diagonal elements correspond to total group areas.
Options:
    -sel <string> [<string>]
        Selection texts for one or more groups.
        Only coordinate-independent selections are supported.
        Selections should not overlap.
)";
    }
protected:

    void before_spawn() override {        
        // Get selection texts
        vector<string> sel_texts = options("sel").as_strings();
        if(sel_texts.size()<1) throw PterosError("At least one selection needed!");
        selections.reserve(sel_texts.size());
        for(const auto&s: sel_texts) selections.emplace_back(system,s);
        // Sanity check
        for(const auto&s: selections){
            if(s.coord_dependent()) throw PterosError("All selections should be coordinate independent!");
        }
        if(check_selection_overlap(selections)) throw PterosError("Selections should not overlap!");
    }

    void pre_process() override {
        // On each instance create voronoi processing class
        voro.create(selections);
    }

    void process_frame(const FrameInfo &info) override {
        voro.compute();
    }

    void post_process(const FrameInfo& info) override {
    }


    void collect_data(const std::vector<std::shared_ptr<TaskBase>>& tasks, int n_frames) override
    {
        // Collect from instances. Put everything into the voro of the first instance
        auto h0 = dynamic_cast<voronoi_par*>(tasks[0].get())->voro;
        for(int i=1; i<tasks.size(); ++i){
            auto h = dynamic_cast<voronoi_par*>(tasks[i].get())->voro;
            for(int g=0; g<h0.num_groups(); ++g){
                h0.get_group(g).total_area += h.get_group(g).total_area;
                h0.get_group(g).total_volume += h.get_group(g).total_volume;
            }
            h0.get_interface_areas() += h.get_interface_areas();
        }

        // Average
        for(int g=0; g<h0.num_groups(); ++g){
            h0.get_group(g).total_area /= double(n_frames);
            h0.get_group(g).total_volume /= double(n_frames);
        }
        h0.get_interface_areas() /= double(n_frames);

        // Output
        h0.write_stats(fmt::format("voronoi_{}",get_id()));
    }

private:
    vector<Selection> selections;
    Voronoi3D voro;
};

CREATE_COMPILED_PLUGIN(voronoi_par)




