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
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <fstream>
#include <map>
#include <set>
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

struct Contact {
    vector<Vector2i> life_fr; // Intervals of frames[first:last]
    vector<Vector2f> life_time; // Intervals of life time [first:last]
    Vector2f energy;
    int num_energy;
    float mean_life_time;
    int num_formed;
};


struct comparator {
    bool operator()(const Vector2i& a, const Vector2i& b) const {
        return (a(0)<b(0)) || (a(0)==b(0) && a(1)<b(1));
    }
};


TASK_SERIAL(contacts)
public:


    string help() override {
        return
R"(Purpose:
    Analyzes contacts between two selections.
Output:

Options:
    -sel1, -sel2
        two selections to compute contacts between.
        Selection should not overlap. An error is thrown is they are.
    -periodic <true|false>, default: false
        Account for periodicity when computing contacts.
    -cutoff, default: 0
        Distance cutoff in nm.
        If cutoff==-1 uses largest sum of VDW radii as a cutoff.
        If cutoff==0  uses the min of VdW and Coulommb cutoffs in the force field (if available).
    -padding, default: 0.1
        Padding added to cutoff in the case of VdW radii (cutoff=-1).
    -transient <true|false>, default: false
        If true the contacts with last for single frame only are recorded.
)";
    }

protected:

    void pre_process() override {
        sel1.modify(system, options("sel1").as_string());
        sel2.modify(system, options("sel2").as_string());
        all.modify(system, "all");

        // Check if selection overlap        
        if(check_selection_overlap({sel1,sel2})) throw PterosError("Selections could not overlap!");

        // Set periodicity
        periodic = options("periodic","false").as_bool();

        /*
        // Removing jumps for both selections if requested
        float unwrap_d = options("unwrap","-1").as_float();
        if(unwrap_d>=0){
            // Add our selections to nojump list
            jump_remover.add_atoms(sel1);
            jump_remover.add_atoms(sel2);
            jump_remover.set_unwrap_dist(unwrap_d);
        }
        */

        // Contacts cutoff
        cutoff = options("cutoff","0").as_float();
        // If zero cutoff given search for maximal sum of VDW distances
        if(cutoff==-1){
            log->info("Computing largest sum VDW distances...");
            float maxd = 0.0, vdw1, vdw2;
            int i,j;
            for(i=0; i<sel1.size(); ++i){                
                vdw1 = sel1.vdw(i);
                for(j=0;j<sel2.size();++j){
                    vdw2 = sel2.vdw(j);
                    if(vdw1+vdw2 > maxd) maxd = vdw1+vdw2;
                }
            }
            log->info("\tLargest sum of VDW distances is {}",maxd);
            // Get padding if given. Default is 0.1
            float pad = options("padding","0.1").as_float();
            cutoff = maxd + pad;
        } else if(cutoff==0) {
            cutoff = system.get_force_field().get_cutoff();
        } else {
            float d = system.get_force_field().get_cutoff();
            if(cutoff!=d) log->warn("Requested cutoff {} is different from cutoff in ff {}!",cutoff,d);
        }

        // Keep transient contacts lasting only 1 frame?
        keep_transient = options("transient","false").as_bool();

        en_f.open(options("en_file",fmt::format("energy_{}.dat",get_id())).as_string());
    }     

    void process_frame(const pteros::FrameInfo &info) override {
        // Search for contacts
        vector<Vector2i> bon;
        vector<float> dist_vec;
        vector<Vector2f> pair_en;

        sel1.apply();
        sel2.apply();

        Vector3i pbc = periodic ? fullPBC : noPBC;
        search_contacts(cutoff,sel1,sel2,bon,dist_vec,true,pbc); // global indexes returned!

        Vector2f total_en(0,0);
        pair_en.resize(bon.size());

        // Get energies if possible
        if(system.force_field_ready()){
            total_en = get_energy_for_list(bon,dist_vec,system, &pair_en);
        }

        // Analyze contacts        
        for(int i=0;i<bon.size();++i){
            // Get contact
            auto& c = bon[i];

            // Sort pair
            sort(c.data(), c.data()+c.size());

            for(int n=0; n<2; ++n){
                // ATOM MAP
                int a = c[n];
                if(atom_map.count(a)){
                    // already present
                    atom_map[a] += 1;
                } else {
                    // new atom
                    atom_map[a] = 1;
                }

                // RES MAP
                int r = all.resindex(c[n]);
                if(res_map.count(r)){
                    // already present
                    res_map[r] += 1;
                } else {
                    // new residue
                    res_map[r] = 1;
                }
            }

            // ATOM CONTACTS
            // See if this contact was already found before
            if(atom_contacts.count(c)){
                // Already present
                auto& cur = atom_contacts[c];
                cur.energy += pair_en[i];
                ++cur.num_energy;
                if( cur.life_fr.back()[1] == info.valid_frame-1 ){
                    // Interval continues
                    cur.life_fr.back()[1] = info.valid_frame;
                    cur.life_time.back()[1] = info.absolute_time;
                } else {
                    // New interval starts
                    cur.life_fr.emplace_back(info.valid_frame,info.valid_frame);
                    cur.life_time.emplace_back(info.absolute_time,info.absolute_time);
                }
            } else {
                // New contact
                Contact cur;
                cur.life_time.emplace_back(info.absolute_time,info.absolute_time);
                cur.life_fr.emplace_back(info.valid_frame,info.valid_frame);
                cur.energy = pair_en[i];;
                cur.num_energy = 1;
                atom_contacts[c] = cur;
            }

            // RESIDUE CONTACTS

            // Get pair of residues
            Vector2i r(all.resindex(c(0)), all.resindex(c(1)));

            //See if this contact was already found before
            if(res_contacts.count(r)){
                // Already present
                auto& cur = res_contacts[r];
                cur.energy += pair_en[i];
                ++cur.num_energy;
                if( cur.life_fr.back()[1] == info.valid_frame-1 ){
                    // Interval continues
                    cur.life_fr.back()[1] = info.valid_frame;
                    cur.life_time.back()[1] = info.absolute_time;
                } else if(info.valid_frame - cur.life_fr.back()[1] > 1) {
                    // New interval starts
                    cur.life_fr.emplace_back(info.valid_frame,info.valid_frame);
                    cur.life_time.emplace_back(info.absolute_time,info.absolute_time);
                }
            } else {
                // New contact
                Contact cur;
                cur.life_time.emplace_back(info.absolute_time,info.absolute_time);
                cur.life_fr.emplace_back(info.valid_frame,info.valid_frame);
                cur.energy = pair_en[i];
                cur.num_energy = 1;
                res_contacts[r] = cur;
            }


        }        

        en_f << info.absolute_time << " " << total_en.sum() << " " << total_en.transpose() << endl;
    }

    void post_process(const pteros::FrameInfo &info) override {
        en_f.close();

        // Analyze atom contacts
        auto it = atom_contacts.begin();
        while (it != atom_contacts.end()) {
            auto& c = it->second;

            if(c.num_energy) c.energy /= (float)c.num_energy;

            c.mean_life_time = 0;
            c.num_formed = 0;

            for(auto& t: c.life_time){
                c.mean_life_time += t(1)-t(0);
                if(t(0)!=t(1) || keep_transient) ++c.num_formed;
            }

            if(c.num_formed == 0 && !keep_transient) {
                it = atom_contacts.erase(it);
            } else {
                c.mean_life_time /= c.num_formed;

                it++;
            }
        }

        // Analyze residue contacts
        it = res_contacts.begin();
        while (it != res_contacts.end()) {
            auto& c = it->second;

            if(c.num_energy) c.energy /= (float)c.num_energy;

            c.mean_life_time = 0;
            c.num_formed = 0;

            for(auto& t: c.life_time){
                c.mean_life_time += t(1)-t(0);
                if(t(0)!=t(1) || keep_transient) ++c.num_formed;
            }

            if(c.num_formed == 0 && !keep_transient) {
                it = atom_contacts.erase(it);
            } else {
                c.mean_life_time /= c.num_formed;

                it++;
            }
        }

        // Output atom contacts
        ofstream f(options("oa",fmt::format("atom_contacts_stats_{}.dat",get_id())).as_string());
        f << "# ATOMS" << endl;
        f << "#i\tj\tn_formed\tlife_t\ten" << endl;
        for(const auto& it: atom_contacts){
            int i = it.first(0);
            int j = it.first(1);
            f << sel1.index(i)+1 << ":" << all.name(i) << ":" << all.resname(i) << "\t"
              << sel2.index(j)+1 << ":" << all.name(j) << ":" << all.resname(j) << "\t"
              << it.second.num_formed << "\t"
              << it.second.mean_life_time << "\t"
              << it.second.energy.sum() << endl;
        }
        f.close();

        // Output residue contacts
        f.open(options("or",fmt::format("res_contacts_stats_{}.dat",get_id())).as_string());
        f << "# RESIDUES" << endl;
        f << "#i\tj\tn_formed\tlife_t\ten" << endl;
        for(const auto& it: res_contacts){
            int i = it.first(0);
            int j = it.first(1);
            f << i << "\t" << j << "\t"
              << it.second.num_formed << "\t"
              << it.second.mean_life_time << "\t"
              << it.second.energy.sum() << endl;
        }
        f.close();

        // Output per atom life time map as pdb file
        float maxv = 0.0;
        for(const auto& it: atom_map) if(maxv<it.second) maxv=it.second;

        for(const auto& it: atom_map){
            all.beta(it.first) = 100.0*it.second/float(info.valid_frame)/maxv;
        }
        all.write(fmt::format("atom_map_{}.pdb",get_id()));

        // Output per residue life time map as pdb file
        maxv = 0.0;
        for(const auto& it: res_map) if(maxv<it.second) maxv=it.second;

        for(const auto& it: res_map){
            all("resindex "+to_string(it.first)).set_beta( 100.0*it.second/float(info.valid_frame)/maxv );
        }
        all.write(fmt::format("res_map_{}.pdb",get_id()));
    }

private:    
    Selection sel1, sel2, all;
    float cutoff;
    bool periodic;
    map<Vector2i,Contact,comparator> atom_contacts;
    map<Vector2i,Contact,comparator> res_contacts;
    bool keep_transient;

    ofstream en_f;

    // Per atom life time heat map
    map<int,float> atom_map;
    // Per residue life time heat map
    map<int,float> res_map;
};


CREATE_COMPILED_PLUGIN(contacts)




