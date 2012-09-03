
/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
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

#ifndef RMSF_H
#define RMSF_H

#include "pteros/analysis/trajectory_processor.h"
//#include "pteros/simulation/simulation.h"

namespace pteros {

struct rmsf_group {
    Selection sel; // running selection

    // Selections for all contained residues. They are not stored in any particular order!
    std::vector<Selection> residue_sel;

    // Selections not containing given residue. Only used for energy calculations
    //vector<Selection> residue_sel_rest;

    // Maps index of residue_sel to rmsf of atoms in this residue
    // rmsf[w][r](j) is rmsf of atom j in residue r for window w
    std::vector< std::map<int,Eigen::VectorXd> > rmsf;
    // Per residue rmsf
    // per_res_rmsf[w][r] is mean rmsf of residue r in window w
    std::vector< std::map<int,double> > per_res_rmsf;

    // Per residue energy
    //std::vector<Eigen::VectorXf> per_res_energy;

    string name;

    // RMSD
    std::vector<double> rmsd;

    // Per group rmsf
    std::vector<double> group_rmsf;
    // Variance per residue
    //Eigen::VectorXf variance;

    // Adds new window
    void add_window();
};

class RMSF: public Consumer
{
public:
    RMSF(Trajectory_processor& proc, Options_tree& opt): Consumer(&proc){
        set_options(opt);
    }

    void set_options(Options_tree& opt){
        options = &opt;
    }

    static void print_help();

    virtual ~RMSF();
    void save_results();

private:
    //------------------------------------
    /// Slots, which are going to be connected with the signals
    /// of trajectory processor
    virtual bool process_frame(const Frame_info& info);
    virtual void pre_process();
    virtual void post_process(const Frame_info& info);
    virtual void window_started(const Frame_info& info);
    virtual void window_finished(const Frame_info& info);
    //------------------------------------
    /// Selections to work with
    std::vector<rmsf_group> groups;

    Selection all;

    /*
        /// Optional topology file for energy calculations
        std::string top_file;
        /// Simulation engine and grid searcher
        Simulation simulation;
        Grid_searcher searcher;
        */

    bool do_rmsd;    
    Options_tree* options;    
};

}

#endif // RMSF_H
