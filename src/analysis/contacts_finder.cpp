/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2014, Semen Yesylevskyy
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

#include "pteros/analysis/contacts_finder.h"
#include "boost/bind.hpp"
#include <sstream>
#include "boost/format.hpp"
#include "boost/foreach.hpp"

#include "json_spirit/json_spirit_reader_template.h"
#include "json_spirit/json_spirit_writer_template.h"


using namespace pteros;
using namespace std;
using namespace Eigen;

void Contacts_finder::create(Trajectory_processor& proc, Options_tree& opt){    
    options = &opt;    

    // Set method and cut-off
    string s;
    if( options->count_options("method") ){
        options->get_option("method") >> s >> dist;
        if(s=="cut_off"){
            method = CUT_OFF;
        } else if(s=="vdw"){
            method = VDW_RADII;
        }
    } else {
        method = CUT_OFF;
        dist = 0.25;
    }
    if(method==VDW_RADII) vdw_gap = dist;

    // Set periodic
    is_periodic = options->get_value<bool>("periodic",false);
}

Contacts_finder::Contacts_finder(Trajectory_processor& proc, Options_tree& opt):
    Consumer(&proc)
{
    create(proc,opt);
}

void Contacts_finder::pre_process(){
    // Set all selection
    all.modify(system,"all");

    // Set selections
    sel_pairs.clear();
    string sel1, sel2;
    BOOST_FOREACH(Options_tree* o, options->get_options("selections")){
        o->get_option("") >> sel1 >> sel2;
        //sel1 = o->get_value<string>("");
        Selections_pair aux;
        sel_pairs.push_back(aux);
        sel_pairs.back().sel1.modify(system,sel1);
        sel_pairs.back().sel2.modify(system,sel2);
    }

    real_time.clear();

    // Prepare staticstics
    for(int i=0; i<sel_pairs.size(); ++i){
        sel_pairs[i].atom_contacts.clear();
        sel_pairs[i].res_contacts.clear();
        sel_pairs[i].per_atom_stats.clear();
        sel_pairs[i].per_res_stats.clear();
        sel_pairs[i].atom_contacts_num.clear();
        sel_pairs[i].res_contacts_num.clear();
        sel_pairs[i].mean_energy = 0.0;
    }

    // If VDW-raddi are used, do some preparation.
    if(method==VDW_RADII){
        VDW.resize(all.size());
        // Assign VDW radii based on element.
        // Radii used in Chimera
        // (http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/vdwtables.html#allatom)
        // are used
        for(int i=0;i<all.size();++i)
            switch(all.Name(i)[0]){
                case 'H': VDW[i] = 0.1; break;
                case 'C': VDW[i] = 0.17; break;
                case 'N': VDW[i] = 0.1625; break;
                case 'O': VDW[i] = 0.149; break; //mean value used
                case 'S': VDW[i] = 0.1782; break;
                case 'P': VDW[i] = 0.1871; break;
                default: VDW[i] = 0.17;
            }

        // Now we'll set the dist for grid searching.
        // Find largest possible distance for contacting atoms
        cout << "Computing largest VDW distance for selection pairs..." << endl;
        float maxd = 0.0;
        int i,j,k;
        for(i=0; i<sel_pairs.size(); ++i)
            for(j=0;j<sel_pairs[i].sel1.size();++j)
                for(k=0;k<sel_pairs[i].sel2.size();++k)
                    if(VDW[sel_pairs[i].sel1.Index(j)]+VDW[sel_pairs[i].sel2.Index(k)] > maxd)
                        maxd = VDW[sel_pairs[i].sel1.Index(j)]+VDW[sel_pairs[i].sel2.Index(k)];
        dist = maxd + vdw_gap;
        cout << "\tLargest distance is " << maxd << " search cut-off is " << dist << endl;
    }
}

void Contacts_finder::process_frame(const Frame_info &info){

    int i,j,n;
    Index_pair at_pair,res_pair;
    Atom_contact_info at_info;
    Residue_contact_info res_info;

    ostringstream ss;

    // Iterators
    Atom_contacts_t::iterator it_at;
    Res_contacts_t::iterator it_res;

    // set of pairs already incremeted on this iteration
    set<Index_pair> is_res_incremented;
    set<int> is_atom1_incremented, is_atom2_incremented;
    set<int> is_res1_incremented, is_res2_incremented;

    bool is_energy = true;

    if(system.force_field_ready())
        is_energy = true;
    else
        is_energy = false;

    // Cycle over all selection pairs
    for(i=0; i<sel_pairs.size(); ++i){
        // Clear list of raw contacts
        clist.clear();
        // Set frame for selections
        sel_pairs[i].sel1.set_frame(0);
        sel_pairs[i].sel2.set_frame(0);
        all.set_frame(0);

        // Search for contacts
        if(method==VDW_RADII){
            // If VDW radii are used we should filter the list first
            // We make aux list and then filter it into clist
            vector<Vector2i> aux;
            Grid_searcher(dist, sel_pairs[i].sel1, sel_pairs[i].sel2, aux, true, is_periodic);
            clist.reserve(aux.size());
            for(int k=0;k<aux.size();++k){
                // See if this contact falls into VDW1+VDW2+gap. If so add to clist
                if( VDW[aux[k](0)]+VDW[aux[k](1)]+vdw_gap >=
                    (all.XYZ(aux[k](0))-all.XYZ(aux[k](1))).norm()
                  ){
                      clist.push_back(aux[k]);
                  }
            }
        } else {
            // Otherwise just call searcher and form clist directly
            //searcher.search(sel_pairs[i].sel1, sel_pairs[i].sel2, clist);
            Grid_searcher(dist, sel_pairs[i].sel1, sel_pairs[i].sel2, clist, true, is_periodic);
        }

        // Add time to the list
        real_time[info.valid_frame] = info.absolute_time;

        // Start parsing contacts list
        n = clist.size();
        is_res_incremented.clear();
        is_atom1_incremented.clear();
        is_atom2_incremented.clear();
        is_res1_incremented.clear();
        is_res2_incremented.clear();

        // Add number of atom contacts for this frame
        sel_pairs[i].atom_contacts_num.push_back(n);

        float total_energy = 0.0;

        // Cycle over all pairs in clist
        for(j=0;j<n;++j){
            //============================
            // Make atom pair
            at_pair.first = clist[j](0);
            at_pair.second = clist[j](1);
            // Compute energy of this pair if requested
            float at_pair_energy;
            if(is_energy){
                //Frame* fp = &system.Frame_data(sel_pairs[0].sel1.get_frame());
                //at_pair_energy = simulation->non_bond_energy(at_pair.first,at_pair.second,*fp).total;
                Energy_components e;
                system.add_non_bond_energy(e,at_pair.first, at_pair.second,
                                           sel_pairs[0].sel1.get_frame(),true);
                at_pair_energy = e.total;
            } else
                at_pair_energy = 0.0;
            total_energy += at_pair_energy;

            // Try to find it in the list of pairs
            it_at = sel_pairs[i].atom_contacts.find(at_pair);
            if(it_at==sel_pairs[i].atom_contacts.end()){
                // No such pair yet, insert new element
                it_at = sel_pairs[i].atom_contacts.insert(make_pair(at_pair,at_info)).first;
                it_at->second.prob = 1;
                it_at->second.mean_energy = 0.0;
                // Add new lifetime entry. It lasts for 1 frame now.
                it_at->second.life_times.push_back(Index_pair(info.valid_frame,1));
                it_at->second.mean_energy = 0.0;
            } else {
                // Pair is already present
                it_at->second.prob += 1;
                // Deal with lifetime.
                // If contact was there last frame, then increase duration
                if(it_at->second.life_times.back().first + it_at->second.life_times.back().second == info.valid_frame)
                    it_at->second.life_times.back().second++;
                else // Otherwise start new interval
                    it_at->second.life_times.push_back(Index_pair(info.valid_frame,1));
            }
            // Add energy to mean energy of this pair
            it_at->second.mean_energy += at_pair_energy;

            //====================================================
            // For both interacting atoms collect some staticstics
            if(sel_pairs[i].per_atom_stats.count(at_pair.first)==0){
                sel_pairs[i].per_atom_stats[at_pair.first].prob = 1;
                is_atom1_incremented.insert(at_pair.first);
            } else {
                // Increment counter only once
                if(is_atom1_incremented.count(at_pair.first)==0){
                    sel_pairs[i].per_atom_stats[at_pair.first].prob += 1;
                    is_atom1_incremented.insert(at_pair.first);
                }
            }
            sel_pairs[i].per_atom_stats[at_pair.first].partners.insert(at_pair.second); // no duplicates

            if(sel_pairs[i].per_atom_stats.count(at_pair.second)==0){
                sel_pairs[i].per_atom_stats[at_pair.second].prob = 1;
                is_atom2_incremented.insert(at_pair.second);
            } else {
                // Increment counter only once
                if(is_atom2_incremented.count(at_pair.second)==0){
                    sel_pairs[i].per_atom_stats[at_pair.second].prob += 1;
                    is_atom2_incremented.insert(at_pair.second);
                }
            }
            sel_pairs[i].per_atom_stats[at_pair.second].partners.insert(at_pair.first); // no duplicates

            //============================
            // Now work with residue pairs

            res_pair = Index_pair(all.Resindex(at_pair.first),all.Resindex(at_pair.second));
            // See if such residue pair exists
            it_res = sel_pairs[i].res_contacts.find(res_pair);
            if(it_res==sel_pairs[i].res_contacts.end()){ // No such pair
                it_res = sel_pairs[i].res_contacts.insert(make_pair(res_pair,res_info)).first;
                it_res->second.prob = 1;
                is_res_incremented.insert(res_pair); // mark as incremented
                it_res->second.mean_energy = 0.0;
            } else { // Pair already exists
                if(is_res_incremented.count(res_pair)==0){
                    it_res->second.prob += 1;
                    is_res_incremented.insert(res_pair); // mark as incremented
                }
            }
            it_res->second.mean_energy += at_pair_energy;

            // Add atom pair to this res pair
            it_res->second.atom_contacts.insert(at_pair);

            //====================================================
            // For both interacting residues collect some staticstics
            if(sel_pairs[i].per_res_stats.count(res_pair.first)==0){
                sel_pairs[i].per_res_stats[res_pair.first].prob = 1;
                is_res1_incremented.insert(res_pair.first);
            } else {
                // Increment counter only once
                if(is_res1_incremented.count(res_pair.first)==0){
                    sel_pairs[i].per_res_stats[res_pair.first].prob += 1;
                    is_res1_incremented.insert(res_pair.first);
                }
            }
            sel_pairs[i].per_res_stats[res_pair.first].partners.insert(res_pair.second);

            if(sel_pairs[i].per_res_stats.count(res_pair.second)==0){
                sel_pairs[i].per_res_stats[res_pair.second].prob = 1;
                is_res2_incremented.insert(res_pair.second);
            } else {
                // Increment counter only once
                if(is_res2_incremented.count(res_pair.second)==0){
                    sel_pairs[i].per_res_stats[res_pair.second].prob += 1;
                    is_res2_incremented.insert(res_pair.second);
                }
            }
            sel_pairs[i].per_res_stats[res_pair.second].partners.insert(res_pair.first);

        }
        // Finished with contacts list

        // Add number of residue contacts for this frame
        sel_pairs[i].res_contacts_num.push_back(is_res_incremented.size());

        // Compute energies
        sel_pairs[i].energy.push_back(total_energy);
        sel_pairs[i].mean_energy += total_energy;//e.total;

    }    
}

void Contacts_finder::post_process(const Frame_info &info){

    // Iterators
    Atom_contacts_t::iterator it;
    Res_contacts_t::iterator ir;

    // Cycle over all selection pairs
    for(int i=0; i<sel_pairs.size(); ++i){

        // For all atom pairs compute probabilities and mean lifetimes
        // Also collect statistics per each interacting atom
        for(it=sel_pairs[i].atom_contacts.begin();it!=sel_pairs[i].atom_contacts.end();it++){
            it->second.prob /= (float)info.valid_frame;
            it->second.mean_energy /= (float)info.valid_frame;

            // Compute mean lifetime
            it->second.mean_lifetime = 0.0;
            for(int j=0;j<it->second.life_times.size();++j)
                it->second.mean_lifetime += it->second.life_times[j].second;
            it->second.mean_lifetime /= (float)it->second.life_times.size();
        }

        // Residue contacts
        for(ir=sel_pairs[i].res_contacts.begin();ir!=sel_pairs[i].res_contacts.end();ir++){
            ir->second.prob /= (float)info.valid_frame;
            ir->second.mean_energy /= (float)info.valid_frame;
        }

        Per_atom_t::iterator ait;
        Per_res_t::iterator rit;
        // Process statistics per each interacting atom and residue
        for(ait=sel_pairs[i].per_atom_stats.begin();ait!=sel_pairs[i].per_atom_stats.end();ait++)
            ait->second.prob /= (float)info.valid_frame;

        for(rit=sel_pairs[i].per_res_stats.begin();rit!=sel_pairs[i].per_res_stats.end();rit++)
            rit->second.prob /= (float)info.valid_frame;

        // Set mean energy
        //if(is_energy){
            sel_pairs[i].mean_energy /= (float)info.valid_frame;
        //}

    } // Sel pair
}

/*
void Contacts_finder::print_info(ostream& out){
    // Iterators
    Atom_contacts_t::iterator it;
    Res_contacts_t::iterator ir;

    int a1,a2;
    string r1,r2;
    float dt = real_time[2]-real_time[1];
    std::set<Index_pair>::iterator cit; // contained atom

    out << "|==============|" << endl;
    out << "| Results:     |" << endl;
    out << "|==============|" << endl;

    // Cycle over all selection pairs
    for(int i=0; i<sel_pairs.size(); ++i){
        out << "-------------------------" << endl;
        out << "  Pair of selections #" << i << endl;
        out << "-------------------------" << endl;
        out << "Atom-Atom contacts:" << endl;
        out << "At1\t\t\tAt2\t\t\tProb\t<life>\t\t#formation" << endl;
        out << "--------------------------------------------------------------------" << endl;
        out << endl;

        for(it=sel_pairs[i].atom_contacts.begin();it!=sel_pairs[i].atom_contacts.end();it++){
            // Dump atom contacts
            a1 = it->first.first;
            a2 = it->first.second;
            out << boost::format("%d.%d.%c.%s.%d\t%d.%d.%c.%s.%d\t%-.3f\t%-.3f (%-.3f)\t%d (%-.3f)")
                    % a1
                    % all.Name(a1)
                    % all.Chain(a1)
                    % all.Resname(a1)
                    % all.Resid(a1)
                    % a2
                    % all.Name(a2)
                    % all.Chain(a2)
                    % all.Resname(a2)
                    % all.Resid(a2)
                    % it->second.prob
                    % it->second.mean_lifetime
                    % (it->second.mean_lifetime*dt)
                    % it->second.life_times.size()
                    % (it->second.life_times.size()*dt)
                << endl;

            // Dump life times
            out << "\tLife times (fr):\t";
            for(int j=0;j<it->second.life_times.size();++j){
                out << boost::format("%d:%d")
                        % it->second.life_times[j].first
                        % (it->second.life_times[j].first + it->second.life_times[j].second)
                    << " ";
            }
            out << endl;
            out << "\tLife times (t):\t\t";
            for(int j=0;j<it->second.life_times.size();++j){
                out << boost::format("%-.3f:%-.3f")
                        % (it->second.life_times[j].first*dt)
                        % ((it->second.life_times[j].first + it->second.life_times[j].second)*dt)
                    << " ";
            }
            out << endl;
        }

        out << "-------------------------" << endl;
        out << "Residue-Residue contacts:" << endl;
        out << "R1\t\t\tR2\t\t\tProb" << endl;
        out << "--------------------------------------------------------------------" << endl;
        // Dump residues
        for(ir=sel_pairs[i].res_contacts.begin();ir!=sel_pairs[i].res_contacts.end();ir++){
            // Dump res contacts
            r1 = ir->first.first;
            r2 = ir->first.second;
            out << boost::format("%s\t%s\t%-.3f") % r1 % r2 % ir->second.prob << endl;
            out << "Contains:" << endl;
            for(cit=ir->second.atom_contacts.begin();cit!=ir->second.atom_contacts.end();cit++){
                a1 = (*cit).first;
                a2 = (*cit).second;
                out << "\t";
                out << boost::format("%d.%d.%c.%s.%d\t%d.%d.%c.%s.%d\t%-.3f")
                    % a1
                    % all.Name(a1)
                    % all.Chain(a1)
                    % all.Resname(a1)
                    % all.Resid(a1)
                    % a2
                    % all.Name(a2)
                    % all.Chain(a2)
                    % all.Resname(a2)
                    % all.Resid(a2)
                    % sel_pairs[i].atom_contacts[*cit].prob
                << endl;
            }
        }
    }
}
*/
// Save whole working session in the json format
void Contacts_finder::save(ostream& out, bool human_readable){

    using namespace json_spirit;

    // Iterators
    Atom_contacts_t::iterator it;
    Res_contacts_t::iterator ir;
    Per_atom_t::iterator ait;
    Per_res_t::iterator rit;

    int a1,a2;
    int r1,r2;
    float dt = real_time[2]-real_time[1];
    std::set<Index_pair>::iterator cit; // contained atom

    mObject result;

    // Cycle over all selection pairs
    mArray selPairs;
    for(int i=0; i<sel_pairs.size(); ++i){
        mObject selPair;
        // Add selections
        selPair["sel1"] = sel_pairs[i].sel1.get_text();
        selPair["sel2"] = sel_pairs[i].sel2.get_text();
        // Add atom contacts
        mArray atomContacts;
        for(it=sel_pairs[i].atom_contacts.begin();it!=sel_pairs[i].atom_contacts.end();it++){
            mObject contact;

            contact["atom1"] = it->first.first;
            contact["atom2"] = it->first.second;
            contact["probability"] = it->second.prob;
            contact["energy"] = it->second.mean_energy;
            contact["mean_lifetime"] = it->second.mean_lifetime;
            contact["times_formed"] = (int)it->second.life_times.size();

            mArray lifeTimes;
            for(int j=0;j<it->second.life_times.size();++j){
                mArray lifeTime;
                lifeTime.push_back(it->second.life_times[j].first);
                lifeTime.push_back(it->second.life_times[j].first + it->second.life_times[j].second);
                lifeTimes.push_back(lifeTime);
            }
            contact["life_times"] = lifeTimes;

            atomContacts.push_back(contact);
        }
        selPair["atom_contacts"] = atomContacts;

        // Residue contacts
        mArray resContacts;
        for(ir=sel_pairs[i].res_contacts.begin();ir!=sel_pairs[i].res_contacts.end();ir++){
            mObject contact;

            contact["residue1"] = ir->first.first;
            contact["residue2"] = ir->first.second;
            contact["probability"] = ir->second.prob;
            contact["energy"] = ir->second.mean_energy;

            // atom contacts
            mArray atomCont;
            mArray arr;
            for(cit=ir->second.atom_contacts.begin();cit!=ir->second.atom_contacts.end();cit++){
                arr.clear();
                arr.push_back((*cit).first);
                arr.push_back((*cit).second);
                atomCont.push_back(arr);
            }
            contact["atom_contacts"] = atomCont;
            resContacts.push_back(contact);
        }
        selPair["residue_contacts"] = resContacts;

        // Per-atom stats
        mArray atomStats;
        for(ait=sel_pairs[i].per_atom_stats.begin();ait!=sel_pairs[i].per_atom_stats.end();ait++){
            mObject atomStat;
            atomStat["atom"] = ait->first;
            atomStat["prob"] = ait->second.prob;

            mArray partn;
            set<int>::iterator t;
            for(t=ait->second.partners.begin();t!=ait->second.partners.end();t++)
                partn.push_back(*t);
            atomStat["partners"] = partn;

            atomStats.push_back(atomStat);
        }
        selPair["per_atom_stats"] = atomStats;

        // Per-residue stats
        mArray resStats;
        for(rit=sel_pairs[i].per_res_stats.begin();rit!=sel_pairs[i].per_res_stats.end();rit++){
            mObject resStat;

            //mObject resRec;
            //resRec["chain"] = rit->first.substr(0,rit->first.find('.'));
            //resRec["resid"] = rit->first.substr(1+rit->first.find('.'));
            resStat["residue"] = rit->first;
            resStat["prob"] = rit->second.prob;

            mArray partn;
            set<int>::iterator t;
            for(t=rit->second.partners.begin();t!=rit->second.partners.end();t++){
                //mObject resRec;
                //resRec["chain"] = t->substr(0,t->find('.'));
                //resRec["resid"] = t->substr(1+t->find('.'));
                partn.push_back(*t);
            }
            resStat["partners"] = partn;

            resStats.push_back(resStat);
        }
        selPair["per_residue_stats"] = resStats;

        // Info about all involved atoms
        mObject usedAtoms;
        for(ait=sel_pairs[i].per_atom_stats.begin();ait!=sel_pairs[i].per_atom_stats.end();ait++){
            mObject at;
            stringstream ss;
            at["name"] = all.Name(ait->first);
            ss << all.Chain(ait->first);
            at["chain"] = ss.str();
            at["resid"] = (int)all.Resid(ait->first);
            at["resname"] = all.Resname(ait->first);

            ss.clear(); ss.str(""); ss << ait->first;
            usedAtoms[ss.str()] = at;
        }
        selPair["atoms_info"] = usedAtoms;

        // Info about all involved residues
        mObject usedRes;
        for(rit=sel_pairs[i].per_res_stats.begin();rit!=sel_pairs[i].per_res_stats.end();rit++){
            // Make selection for first atom in this residue to retrieve info
            /*
            stringstream ss; ss << "chain " << rit->first.substr(0,rit->first.find('.'))
                                << " and "
                                << "resid " << rit->first.substr(1+rit->first.find('.'));

            int ind = Selection(*all.get_system(),ss.str()).get_index()[0];
            */
            stringstream ss;
            ss << "resindex " << rit->first;
            Selection tmp_sel(*all.get_system(),ss.str());
            int ind = tmp_sel.get_index()[0];
            mObject r;
            r["chain"] = all.Chain(ind);
            r["resid"] = (int)all.Resid(ind);
            r["resname"] = all.Resname(ind);

            ss.clear(); ss.str(""); ss << rit->first;
            usedRes[ss.str()] = r;
        }
        selPair["residues_info"] = usedRes;

        selPair["mean_energy"] = sel_pairs[i].mean_energy;

        /*
        // Number of atom contacts with time
        mArray num;
        for(int j=0;j<sel_pairs[i].atom_contacts_num.size();++j)
            num.push_back(sel_pairs[i].atom_contacts_num[j]);
        selPair["atom_contacts_num"] = num;

        // Number of residue contacts with time
        num.clear();
        for(int j=0;j<sel_pairs[i].res_contacts_num.size();++j)
            num.push_back(sel_pairs[i].res_contacts_num[j]);
        selPair["res_contacts_num"] = num;
        */

        selPairs.push_back(selPair);
    }
    result["time_step"] = dt;
    result["time_unit"] = "ps";
    result["selection_pairs"] = selPairs;
    result["is_periodic"] = is_periodic;
    switch(method){
        case CUT_OFF:
            result["method"] = "cut-off";
            result["distance"] = dist;
            break;
        case VDW_RADII:
            result["method"] = "VDW-radii";
            result["distance"] = vdw_gap;
            break;
    }

    result["is_energy"] = true;//(simulation!=NULL);

    // Add info for trajectories
    result["job_info"] = options->to_json_string();

    write_stream(mValue(result),out,human_readable);
}

// Write statistics per frame
void Contacts_finder::per_frame_stats(ostream& out){
    int i;
    out << "# fr time (atom_contacts_num res_contacts_num energy)*n" << endl;
    int num = sel_pairs[0].energy.size();

    for(int fr=0; fr<num; ++fr){
        out << fr << " " << real_time[fr] << " ";

        for(i=0;i<sel_pairs.size();++i){
            out << sel_pairs[i].atom_contacts_num[fr] << " ";
            out << sel_pairs[i].res_contacts_num[fr] << " ";
            out << sel_pairs[i].energy[fr] << " ";
        }
        out << endl;
    }
}

void Contacts_finder::print_help(){
    cout << "Specific options:\n"
            "-----------------------------------------------------\n"
            "--method [cut_off|vdw] distance\n"
            "\tOptional. Method of contacts definition and contact distance in nm.\n"
            "\tIf vdw is specified, then the distance is added to VDW radii of two contacting atoms.\n"
            "\tDefaults to 'cut_off 0.25'\n"

            "--selections sel1 sel2:\n"
            "\tRequired. Pair of selections. The contacts are found between them.\n"
            "\tThis option could be given several times.\n"

            "\n"
         << endl;
}
