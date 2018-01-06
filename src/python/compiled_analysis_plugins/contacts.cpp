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


#include "pteros/python/compiled_plugin.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <fstream>
#include <map>

using namespace std;
using namespace pteros;
using namespace Eigen;

struct Contact {
    vector<Vector2i> life_fr; // Intervals of frames[first:last]
    vector<Vector2f> life_time; // Intervals of life time [first:last]
    float energy;
    int num_energy;
    float mean_life_time;
    int num_formed;
};


struct comparator {
    bool operator()(const Vector2i& a, const Vector2i& b) const {
        return (a(0)<b(0)) || (a(0)==b(0) && a(1)<b(1));
    }
};

class contacts: public pteros::Compiled_plugin_base {
public:
    contacts(pteros::Trajectory_processor* pr, const pteros::Options& opt): Compiled_plugin_base(pr,opt) {}

    string help(){
        return  "Purpose:\n"                
                ;
    }

protected:

    void pre_process(){
        sel1.modify(system, options("sel1").as_string());
        sel2.modify(system, options("sel2").as_string());

        // Check if selection overlap
        {
            vector<int> v;
            set_intersection(sel1.index_begin(),sel1.index_end(),
                             sel2.index_begin(),sel2.index_end(),v.begin());
            if(v.size()) throw Pteros_error("Selections could not overlap!");
        }

        // Removing jumps for both selections if requested
        float unwrap_d = options("nojump","0").as_float();
        if(unwrap_d>=0){
            // Add our selections to nojump list
            jump_remover.add_atoms(sel1);
            jump_remover.add_atoms(sel2);
            jump_remover.set_unwrap_dist(unwrap_d);
        }

        // Contacts cutoff
        cutoff = options("cutoff","0").as_float();
        // If zero cutoff given search for maximal sum of VDW distances
        if(cutoff==0){
            cout << "Computing largest sum VDW distances..." << endl;
            float maxd = 0.0, vdw1, vdw2;
            int i,j;
            for(i=0; i<sel1.size(); ++i){
                vdw1 = sel1.VDW(i);
                for(j=0;j<sel2.size();++j){
                    vdw2 = sel2.VDW(i);
                    if(vdw1+vdw2 > maxd) maxd = vdw1+vdw2;
                }
            }
            cout << "\tLargest sum of VDW distances is " << maxd << endl;
            // Get padding if given. Default is 0.1
            float pad = options("padding","0.1").as_float();
            cutoff = maxd + pad;
        }
        cout << "Search distances is " << cutoff << endl;

        // Set periodicity
        periodic = options("periodic","false").as_bool();

        // Keep transient contacts lasting only 1 frame?
        keep_transient = options("transient","false").as_bool();

        en_f.open(options("en_file","energy_"+label+".dat").as_string());
    }     

    void process_frame(const pteros::Frame_info &info){                
        // Search for contacts
        vector<Vector2i> bon;
        vector<float> dist_vec;
        search_contacts(cutoff,sel1,sel2,bon,false,periodic,&dist_vec); // Local indexes returned
        // Analyze contacts
        float total_en = 0.0;
        for(auto c: bon){
            // Sort pair
            sort(c.data(), c.data()+c.size());

            // ENERGY

            float en = 0.0;

            if(system.force_field_ready()){
                en = system.non_bond_energy(c(0),c(1),0,periodic).total;
            }

            total_en += en;

            // ATOM CONTACTS

            // See if this contact was already found before
            if(atom_contacts.count(c)){
                // Already present
                auto& cur = atom_contacts[c];
                cur.energy += en;
                ++cur.num_energy;
                if( cur.life_fr.back()[1] == info.valid_frame-1 ){
                    // Interval continues
                    cur.life_fr.back()[1] = info.valid_frame;
                    cur.life_time.back()[1] = info.absolute_time;
                } else {
                    // New interval starts
                    cur.life_fr.push_back(Vector2i(info.valid_frame,info.valid_frame));
                    cur.life_time.push_back(Vector2f(info.absolute_time,info.absolute_time));
                }
            } else {
                // New contact
                Contact cur;
                cur.life_time.push_back(Vector2f(info.absolute_time,info.absolute_time));
                cur.life_fr.push_back(Vector2i(info.valid_frame,info.valid_frame));
                cur.energy = en;
                cur.num_energy = 1;
                atom_contacts[c] = cur;
            }

            // RESIDUE CONTACTS

            // Get pair of residues
            Vector2i r(sel1.Resindex(c(0)), sel2.Resindex(c(1)));

            //See if this contact was already found before
            if(res_contacts.count(r)){
                // Already present
                auto& cur = res_contacts[r];
                cur.energy += en;
                ++cur.num_energy;
                if( cur.life_fr.back()[1] == info.valid_frame-1 ){
                    // Interval continues
                    cur.life_fr.back()[1] = info.valid_frame;
                    cur.life_time.back()[1] = info.absolute_time;
                } else if(info.valid_frame - cur.life_fr.back()[1] > 1) {
                    // New interval starts
                    cur.life_fr.push_back(Vector2i(info.valid_frame,info.valid_frame));
                    cur.life_time.push_back(Vector2f(info.absolute_time,info.absolute_time));
                }
            } else {
                // New contact
                Contact cur;
                cur.life_time.push_back(Vector2f(info.absolute_time,info.absolute_time));
                cur.life_fr.push_back(Vector2i(info.valid_frame,info.valid_frame));
                cur.energy = en;
                cur.num_energy = 1;
                res_contacts[r] = cur;
            }


        }

        en_f << info.absolute_time << " " << total_en << endl;
    }

    void post_process(const pteros::Frame_info &info){
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
        ofstream f(options("oa","contacts_"+label+".dat").as_string());
        f << "# ATOMS" << endl;
        f << "#i\tj\tn_formed\tlife_t\ten" << endl;
        for(const auto& it: atom_contacts){
            int i = it.first(0);
            int j = it.first(1);
            f << sel1.Index(i)+1 << ":" << sel1.Name(i) << ":" << sel1.Resname(i) << "\t"
              << sel2.Index(j)+1 << ":" << sel2.Name(j) << ":" << sel2.Resname(j) << "\t"
              << it.second.num_formed << "\t"
              << it.second.mean_life_time << "\t"
              << it.second.energy << endl;
        }
        f.close();

        // Output residue contacts
        f.open(options("or","rescontacts_"+label+".dat").as_string());
        f << "# RESIDUES" << endl;
        f << "#i\tj\tn_formed\tlife_t\ten" << endl;
        for(const auto& it: res_contacts){
            int i = it.first(0);
            int j = it.first(1);
            f << i << "\t" << j << "\t"
              << it.second.num_formed << "\t"
              << it.second.mean_life_time << "\t"
              << it.second.energy << endl;
        }
        f.close();

    }

private:    
    Selection sel1, sel2;
    float cutoff;
    bool periodic;
    map<Vector2i,Contact,comparator> atom_contacts;
    map<Vector2i,Contact,comparator> res_contacts;
    bool keep_transient;

    ofstream en_f;
};


CREATE_COMPILED_PLUGIN(contacts)

