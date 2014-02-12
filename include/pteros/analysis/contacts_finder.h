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

#ifndef CONTACTS_FINDER_H
#define CONTACTS_FINDER_H

#include <iostream>
#include <string>
#include <set>
#include "pteros/pteros.h"
#include "pteros/analysis/trajectory_processor.h"
#include "pteros/analysis/options_parser.h"
#include "pteros/core/grid_search.h"

namespace pteros {

enum Contacts_def {CUT_OFF, VDW_RADII};

typedef std::pair<int,int> Index_pair;

struct Atom_contact_info {
    float prob; // Probability of formation
    std::vector<Index_pair> life_times; // First value - formation, second - duration
    float mean_lifetime;
    float mean_energy;
};

struct Residue_contact_info {
    float prob; // Probability of formation
    std::set<Index_pair> atom_contacts; // Indexes of atom contacts
    float mean_energy;
};

struct Per_atom_stat {
    float prob;
    std::set<int> partners;
};

struct Per_residue_stat {
    float prob;
    std::set<int> partners;
};

typedef std::map<Index_pair,Atom_contact_info> Atom_contacts_t;
typedef std::map<Index_pair,Residue_contact_info> Res_contacts_t;
typedef std::map<int,Per_atom_stat> Per_atom_t;
typedef std::map<int,Per_residue_stat> Per_res_t;

struct Selections_pair {
    Selection sel1;
    Selection sel2;
    Atom_contacts_t atom_contacts;
    Res_contacts_t res_contacts;
    Per_atom_t per_atom_stats;
    Per_res_t per_res_stats;
    vector<int> atom_contacts_num;
    vector<int> res_contacts_num;
    vector<float> energy;
    float mean_energy;
};

class Contacts_finder: public Consumer {
    public:
        Contacts_finder(Trajectory_processor& proc, Options_tree& opt);

        void create(Trajectory_processor& proc, Options_tree& opt);

        /// Prints pretty human-readable summary to the stream
        //void print_info(std::ostream& out);
        /// Saves results in JSON format to the stream
        void save(std::ostream& out, bool human_readable=true);
        /// Writes statisticts per frame (number of contacts, etc.)
        void per_frame_stats(std::ostream& out);
        /// Print help
        static void print_help();

    private:
        Options_tree* options;

        //------------------------------------
        /// Slots, which are going to be connected with the signals
        /// of trajectory processor
        void process_frame(const Frame_info& info);
        void pre_process();
        void post_process(const Frame_info& info);
        //------------------------------------

        // Grid searcher object
        Grid_searcher searcher;
        // Aux list of raw contacts
        std::vector<Eigen::Vector2i> clist;

        // Mapping of real time
        std::map<int,float> real_time;

        // VDW radii
        std::vector<float> VDW;

        // Pairs of selections to process
        std::vector<Selections_pair> sel_pairs;
        // Parameters
        double dist; // Cut-off for grid searching
        double vdw_gap; // Gap for VDE mode
        bool is_periodic;
        Contacts_def method;

        // All selection
        Selection all;
};

}
#endif
