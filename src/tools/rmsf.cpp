/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2020, Semen Yesylevskyy
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




#include "rmsf.h"

using namespace pteros;
using namespace std;

void rmsf_group::add_window(){
    std::map<int,Eigen::VectorXd> v;
    rmsf.push_back(v);
    // Now allocate size for each residue
    for(int i=0; i<residue_sel.size(); ++i){
        int resind = residue_sel[i].Resindex(0);
        rmsf.back()[resind].resize(residue_sel[i].size());
        rmsf.back()[resind].fill(0.0);
    }

    // Allocate per res
    std::map<int,double> a;
    per_res_rmsf.push_back(a);
    for(int i=0; i<residue_sel.size(); ++i){
        int resind = residue_sel[i].Resindex(0);
        per_res_rmsf.back()[resind] = 0.0;
    }

    /*
    // Allocate per res energy
    per_res_energy.push_back(a);
    per_res_energy.back().resize(residue_sel.size());
    per_res_energy.back().fill(0.0);
    */

    // Add group rmsf
    group_rmsf.push_back(0.0);
}


RMSF::~RMSF()
{
    //dtor
}

void RMSF::pre_process(){
    // Add selections
    int k = 0;
    for(auto sel: options("selections").as_strings()){
        rmsf_group gr;        
        gr.sel.modify(system,sel);
        groups.push_back(gr);
        cout << "Added selection '" << sel << "'" << endl;
        k++;
    }

    /*
    // If topology file is provided set it
    if(options.count_options("topology_file")){
        top_file = options.get_value<string>("topology_file");
        cout << "Creating simulation object..." << endl;
        // Do not run setup() - we don't need searcher and nblist stuff
        simulation.load_topology(*system,top_file);
        // Init grid searcher
        float dist = options.get_value<float>("energy_cut_off",1.0);
        bool is_periodic = options.get_value<bool>("is_periodic",false);
        searcher.create(dist,*system,true,is_periodic);
    } else
        top_file = ""; // No topology
    */

    // For each selection create residue selections and start first window
    for(auto& gr: groups){
        gr.rmsf.clear();
        gr.sel.each_residue(gr.residue_sel);

        /*
        cout << gr.sel.get_text() << " contains " << gr.residue_sel.size() << " residues" << endl;
        for(int i=0; i<gr.residue_sel.size();++i){
            cout << "Res." << gr.residue_sel[i].Resindex(0) << ": ";
            for(int j=0; j<gr.residue_sel[i].size();++j)
                cout << gr.residue_sel[i].Index(j) << " ";
            cout << endl;
        }
        */

        /*
        // Create complementary selections if requested
        if(top_file.size()){
            for(int i=0;i<gr.residue_sel.size();++i){
                gr.residue_sel_rest.push_back(
                    Selection(*system,"not ("+gr.residue_sel[i].get_text()+")")
                );
            }
        }
        */

        gr.add_window();
        // Allocate variance
        //gr.variance.resize(gr.residue_sel.size());
    }

    all.modify(system,"all");

    do_rmsd = options("do_rmsd","false").as_bool();
}

void RMSF::window_started(const Frame_info &info){            
    if(system.num_frames()==1) return; // In case this is the very beginning
    system.frame_copy(0,1);
}

void RMSF::window_finished(const Frame_info &info){

    for(auto& gr: groups){

        int resind;
        // Normilize rmsf for each atom
        for(int r=0; r<gr.residue_sel.size(); ++r){
            resind = gr.residue_sel[r].Resindex(0);
            gr.rmsf[info.win_num][resind] = ( gr.rmsf[info.win_num][resind].array()/
                                        (info.win_last_time - info.win_start_time) ).sqrt();
        }

        // Get per-residue and per group averages
        double av = 0.0, ma = 0.0;
        for(int r=0; r<gr.residue_sel.size(); ++r){
            resind = gr.residue_sel[r].Resindex(0);

            vector<float> tmp = gr.residue_sel[r].get_mass();
            Eigen::Map<Eigen::VectorXf> mm(tmp.data(),tmp.size());

            double v = ( gr.rmsf[info.win_num][resind].array() *
                         mm.array().cast<double>() ).sum();
            double m = mm.sum();
            gr.per_res_rmsf[info.win_num][resind] = v/m;
            av += v;
            ma += m;
        }
        gr.group_rmsf[info.win_num] = av/ma;

        // Add new window to rmsf
        gr.add_window();
    }

    cout << "Window " << info.win_num << " done!" << endl;
}


void RMSF::process_frame(const Frame_info &info){
    // Very first frame should be used as first reference
    if(info.last_frame==info.first_frame){
        cout << "Setting first trajectory frame as initial reference..." << endl;
        if(system.num_frames()==1){
            system.frame_dup(0);
        }        
        system.frame_copy(0,1);
    }

    all.set_frame(0);    

    Eigen::Affine3f t;

    // For each group
    for(auto& gr: groups){
        // Update selection
        gr.sel.set_frame(0);

        // Allign whole structure by this selection
        t = gr.sel.fit_transform(0,1);
        all.apply_transform(t);
        //gr.sel.apply_transform(t);

        if(do_rmsd) gr.rmsd.push_back( gr.sel.rmsd(0,1) );

        // Cycle for all residues
        //cout << ">>" << gr.residue_sel.size() << " " << info.win_num << endl;
        for(int r=0; r<gr.residue_sel.size(); ++r){            
            // For all atoms in selection compute difference with reference
            int resind = gr.residue_sel[r].Resindex(0);

            //cout << r << " " << resind << endl;
            // For each atom in residue r (resid)
            for(int i=0; i<gr.residue_sel[r].size(); ++i){
                gr.rmsf[info.win_num][resind](i) += ( gr.residue_sel[r].XYZ(i,0) -
                                            gr.residue_sel[r].XYZ(i,1) ).squaredNorm();
            }
        }

        /*
        // Run energy calculation if requested
        if(top_file.size()){

            // Generate clist for given residue
            vector<Eigen::Vector2i> clist;
            for(int r=0; r<gr.residue_sel_rest[r].size(); ++r){
                searcher.search(gr.residue_sel[r], gr.residue_sel_rest[r], clist);
                Energy_components e = simulation.non_bond_energy(clist);
                gr.per_res_energy[info.win_num](r) += e.total;
                cout << "==> " << r << " " << e.total << endl;
            }
        }
        */

    }    
}

void RMSF::post_process(const Frame_info &info){
}

void RMSF::save_results(){
    // For each group open file and write matrices
    cout << "Saving data..." << endl;

    string prefix = options("output_prefix","").as_string();

    int i = 0;
    for(auto& gr: groups){
        vector<int> resinds = gr.sel.get_unique_resindex();
        int resind;

        // Per-residue matrix
        ofstream f((prefix+gr.name+"_per_residue_matrix.dat").c_str());
        for(int w=0; w<gr.rmsf.size()-1; ++w){ // cycle over windows
            // cycle over residues
            for(int r=0; r<resinds.size(); ++r){
                f << gr.per_res_rmsf[w][resinds[r]] << " ";
            }
            f << endl;
        }
        f.close();

        // Per-group graph
        f.open((prefix+gr.name+"_per_group_rmsf.dat").c_str());
        for(int w=0; w<gr.rmsf.size()-1; ++w){
            f << w << " " << gr.group_rmsf[w] << endl;
        }
        f.close();
        /*
        // Per-residue time variance
        f.open(("per_residue_variance_"+boost::boost::lexical_cast<string>(i)+".dat").c_str());
        for(int j=0; j<gr.residue_sel.size(); ++j){
            f << j << " " << gr.variance(j) << endl;
        }
        f.close();
        */

        // Per-atom matrix
        // Only save it if requested because ti is HUGE!
        if(options("save_per_atom","false").as_bool()){
            f.open((prefix+gr.name+"_per_atom_rmsf.dat").c_str());
            for(int w=0; w<gr.rmsf.size()-1; ++w){
                // For each group do
                // We need to output residues in ascending order, while
                // they are stored as a map with no particular order
                for(int r=0; r<resinds.size(); ++r){
                    f << gr.rmsf[w][resinds[r]].transpose() << " ";
                }
                f << endl;
            }
            f.close();
        }

        if(do_rmsd){
            f.open((prefix+gr.name+"_rmsd.dat").c_str());
            for(int w=0; w< gr.rmsd.size(); ++w){
                f << gr.rmsd[w] << endl;
            }
            f.close();
        }

        ++i;
    }
}

void RMSF::print_help(){
    cout << "Specific options for tRMSF calculation:\n"
            "---------------------------------------\n"

            "--selections sel_text [--name sel_name]:\n"
            "\tRequired. Selection for doing tRMSF.\n"
            "\tsel_name is used for output files. If not specified numbers are used.\n"
            "\tThis option could be given several times.\n"

            "--do_rmsd [true|false]:\n"
            "\tOptional. Also compute RMSD for each selection.\n"
            "\tIf no window is specified for trajectory processing, then gives normal RMSD.\n"
            "\tOtherwise gives RMSD from the first frame in each window.\n"

            "--save_per_atom [true|false]:\n"
            "\tOptional. If true tRMSF matrix for atoms is written in addition.\n"
            "\tThis file may be HUGE! Defaults to false.\n"

            "--output_prefix str:\n"
            "\tOptional prefix for all output files. May contain the common path for example.\n"

            "\n"
         << endl;
}


