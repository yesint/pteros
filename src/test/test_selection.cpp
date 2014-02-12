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

#include "pteros/pteros.h"
#include "pteros/core/grid_search.h"
#include "pteros/analysis/trajectory_processor.h"
#include "pteros/analysis/rmsf.h"
#include "pteros/core/mol_file.h"
#include <string>

#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace pteros;
using namespace Eigen;


int main(int argc, char** argv)
{
    try{
        /*
        System s1("/home/semen/install/dssp-2.1.0/2lao.pdb");
        Selection all(s1,"resname ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR");
        //Selection all(s1,"resname ALA");
        all.write("prot.pdb");
        s1.dssp("dssp_out.dat");
        cout << "Code string:" << endl;
        cout << s1.dssp() << endl;

        return 0;
        */


        //System sys("/media/data/semen/trajectories/asymmetric_bicelle/no_restr/last.pdb");
        System sys("/media/data/semen/trajectories/grand_challenge/dump.gro");

        Selection sel1(sys,"name ROH");
        Selection sel2;
        sel2 = sel1;
        cout << sel1.size() << " " << sel2.size() << endl;
        cout << sel1.get_text() << endl;
        cout << sel2.get_text() << endl;
        //System sys("/home/semen/work/Projects/pteros/pteros-git-build/release/bin/dppc_vesicle.gro");
        //System sys("/home/semen/work/Projects/pteros/pteros-git-build/release/bin/dpc-vesicle.gro");
        //System sys("/home/semen/work/Projects/pteros/pteros-git-build/release/bin/pope_vesicle.gro");

//        Selection(sys,"resname DOPC DOPS").unwrap_bonds(0.6,Vector3i(1,0,1));
        //Selection(sys,"resname DPPC").unwrap_bonds(0.6,Vector3i(1,0,1));
        cout << "Unwrapping done" << endl;
        Selection(sys,"resname DOPC DOPS").write("wrapped.pdb");

/*
        Matrix3f b = sys.Box(0).get_box();
        b.col(0) *= 2.0;
        b.col(2) *= 2.0;
        sys.Box(0).modify(b);

        Selection(sys,"all").set_mass(1); // CG atoms may have no mass assigned

        Lipid_assembly a;
        a.create(sys,"resname DOPC DOPS","name PO4 NC3 CNO", "name C5A C5B",2.0);
        //a.create(sys,"resname DOPC DOPS","name C5A C5B","name PO4",3.5);

*/

// Middle: 139.699997, 3.000000, 64.099998
// Edge    152.800003, 0.200000, 166.800003



    } catch(const Pteros_error& e){ e.print(); }

    return 0;
    //---------------------------------------

    try{
        System t0("/home/semen/work/Projects/pteros/pteros_git_build/experimental/release/bin/topol.tpr.pttop");
        return 1;

        System t("/media/data/semen/trajectories/RC/simulation_C/after_lip_equil.pdb");
        vector<Eigen::Vector2i> bon;
        clock_t t1,t2;
        float cutoff;

        // Test 1
        cout << "within <cutoff> of resname BCL" << endl;
        cutoff = 0.2;
        while(cutoff<=2.1){
            t1 = clock();
            for(int i=0;i<100;++i){
                Selection(t,"within "+boost::lexical_cast<string>(cutoff)+" of resname BCL");
            }
            t2 = clock();
            cout << cutoff << ": " << (float)(t2-t1)/float(CLOCKS_PER_SEC) << endl;
            cutoff += 0.2;
        }

        return 1;

        // Test 2
        cout << "resname SOL and within <cutoff> of resname SOL" << endl;
        cutoff = 0.2;
        while(cutoff<=2.1){
            t1 = clock();
            for(int i=0;i<100;++i){
                Selection(t,"resname SOL and within "+boost::lexical_cast<string>(cutoff)+" of resname SOL");
            }
            t2 = clock();
            cout << cutoff << ": " << (float)(t2-t1)/float(CLOCKS_PER_SEC) << endl;
            cutoff += 0.2;
        }

        // Test 3
        cout << "name OW and within <cutoff> of resname ALA" << endl;
        cutoff = 0.2;
        while(cutoff<=2.1){
            t1 = clock();
            for(int i=0;i<100;++i){
                Selection(t,"name OW and within "+boost::lexical_cast<string>(cutoff)+" of resname ALA");
            }
            t2 = clock();
            cout << cutoff << ": " << (float)(t2-t1)/float(CLOCKS_PER_SEC) << endl;
            cutoff += 0.2;
        }

        // Test 4 - creation of coord-independent selections
        cout << "resid <i> and name CA" << endl;
        cutoff = 0.2;
            t1 = clock();
            for(int i=0;i<500;++i){
                Selection(t,"resid "+boost::lexical_cast<string>(i)+" and name CA");
            }
            t2 = clock();
            cout << ": " << (float)(t2-t1)/float(CLOCKS_PER_SEC) << endl;


        // Test 5 - creation of coord-independent selections by residue
        cout << "by residue resid <i> and name CA" << endl;
        cutoff = 0.2;
            t1 = clock();
            for(int i=0;i<500;++i){
                Selection(t,"by residue resid "+boost::lexical_cast<string>(i)+" and name CA");
            }
            t2 = clock();
            cout << ": " << (float)(t2-t1)/float(CLOCKS_PER_SEC) << endl;



        /*
        ofstream f("out.dat");
        for(int fr=1;fr<t.num_frames();++fr){
            Energy_components e;
            t.add_non_bond_energy(e,0,1,fr,true);
            f << fr << "  " << e.to_str() << endl;
        }
        f.close();
        */

        /*
        Vector3f ref = s.XYZ(10000);
        for(int i=0;i<s.size();++i){
            s.XYZ(i) = t.get_closest_image(s.XYZ(i),ref,0);
        }
        s.write("/home/semen/work/Projects/pteros/pteros_git/src/test/data/2lao-tr.gro");
        */

        /*
        cout << s.get_box_volume() << endl;
        Vector3f moments;
        Matrix3f axes;
        s.inertia(moments,axes);
        cout << moments << endl;
        cout << axes << endl;
        */

        //Affine3f tr = s.principal_transform(false);
        //s.apply_transform(tr);
        //s.write("/home/semen/work/Projects/pteros/pteros_git/src/test/data/2lao-tr.gro");

        //cout << s.center(true,true) << endl;
        //s.wrap();
        //s.unwrap();


/*

        ifstream f("inp.indent");
        string str((istreambuf_iterator<char>(f)), istreambuf_iterator<char>());
        f.close();
        Options_tree opt;
        opt.from_indented(str);
        cout << opt.to_json_string() << endl;
        cout << opt.to_indented() << endl;

        return 0;
*/
        /*
        Options_tree options;
        options.from_command_line(argc,argv);

        Trajectory_processor proc(options);
        boost::shared_ptr<Consumer> c1(new Consumer(&proc));
        boost::shared_ptr<Consumer> c2(new Consumer(&proc));
        boost::shared_ptr<Consumer> c3(new Consumer(&proc));
        vector<string> sv;
        sv.push_back(s);
        proc.compute();
        */
        //string fname("/home/semen/work/Projects/kornelyuk/Sasha/dimer_md/1/dimer_pdb2gmx.gro");
        //string fname("/home/semen/work/Projects/asymmetric_bilayer/prepare_small/with_ions.pdb");


        //Selection sel(sys,"not resname W ION");
        //Bilayer bi(sel,"name PO4",2.0);

        //Selection chol(sys,"resname CHOL");
        //Bilayer_point_info info = bi.point_info(chol.XYZ(0));
        //info.print();



        //Selection sel(sys,"within 0.6 of (resname BCL and name MG)");
        //cout << sel.get_text() << endl;

/*
        Selection sel(sys,"all");
        Selection sel1(sys,"resid 1-300");
        Selection sel2(sys,"resid 302-600");
        vector<Vector2i> bon;



        Grid_searcher s;
        s.assign_to_grid(0.5,sel,true,false);
        vector<int> l;

        clock_t t1 = clock();
        for(int i=0; i<0; ++i){
            //s.search_within(sel1,l);
            Selection sel3(sys,"within  of (resid 1-300)");
            //s.search_within(sel.XYZ(100),l);
            //Grid_searcher(0.26,sel,bon,true);
            //Grid_searcher(0.5,sel1,sel2,bon,true,true);
        }

        cout << "Elapsed:" << (clock()-t1)/float(CLOCKS_PER_SEC) << endl;

        //Selection sel3(sys,"dist pbc point 0.3 0.2 0.7> 3");
        Selection sel3(sys,"dist pbc vector 0.3 0.2 0.7 1 -1.0 1 > 4");
        //Selection sel3(sys,"within 3.5 periodic of (resid 313 and name PO4)");
        sel3.set_chain('Q');
        Selection(sys,"all").write("sel.pdb");



        cout << sel3.size() << endl;
        cout << l.size() << endl;
        //for(int i=0;i<l.size();++i) cout << l[i] << " ";
        cout << endl;
*/
/*
        sel1.X(0) =  -2.803;
        sel1.Y(0) = -15.373;
        sel1.Z(0) =  24.556;
        sel1.X(1) =   0.893;
        sel1.Y(1) = -16.062;
        sel1.Z(1) =  25.147;
        sel1.X(2) =   1.368;
        sel1.Y(2) = -12.371;
        sel1.Z(2) =  25.885;
        sel1.X(3) =  -1.651;
        sel1.Y(3) = -12.153;
        sel1.Z(3) =  28.177;
        sel1.X(4) =  -0.440;
        sel1.Y(4) = -15.218;
        sel1.Z(4) =  30.068;
        sel1.X(5) =   2.551;
        sel1.Y(5) = -13.273;
        sel1.Z(5) =  31.372;
        sel1.X(6) =   0.105;
        sel1.Y(6) = -11.330;
        sel1.Z(6) =  33.567;

        sel2.X(0) = -14.739;
        sel2.Y(0) = -18.673;
        sel2.Z(0) =  15.040;
        sel2.X(1) = -12.473;
        sel2.Y(1) = -15.810;
        sel2.Z(1) =  16.074;
        sel2.X(2) = -14.802;
        sel2.Y(2) = -13.307;
        sel2.Z(2) =  14.408;
        sel2.X(3) = -17.782;
        sel2.Y(3) = -14.852;
        sel2.Z(3) =  16.171;
        sel2.X(4) = -16.124;
        sel2.Y(4) = -14.617;
        sel2.Z(4) =  19.584;
        sel2.X(5) = -15.029;
        sel2.Y(5) = -11.037;
        sel2.Z(5) =  18.902;
        sel2.X(6) = -18.577;
        sel2.Y(6) = -10.001;
        sel2.Z(6) =  17.996;

        sel1.set_mass(1.0);
        sel2.set_mass(1.0);

        /*
        for(int i=0; i<con.size(); ++i){
            cout << con[i](0) << " " << con[i](1) << " "
                 << (sys.XYZ(con[i](1),1)-sys.XYZ(con[i](0),1)).norm() << endl;
        }
        */


        /*
        System s;
        s.load(path+"2lao.pdb");
        Selection sel1(s,"name CA and resid 2  to 5");
        cout << sel1.size() << endl;


        System s;
        s.load(path+"2lao.pdb");
        s.frame_dup(0);
        Selection sel1(s,"name CA");
        Selection sel2(s,"name CA");
        sel2.set_frame(1);
        sel2.rotate(1,0.3);
        cout << rmsd(sel1,0,sel2,1) << endl;
        Affine3f t = fit_transform(sel1,sel2);
        sel1.apply_transform(t);
        cout << rmsd(sel1,0,sel2,1) << endl;
        */

/*
    cout << "==================================="<<endl;
    cout << "Testing constructors of Selection" << endl;
    cout << "==================================="<<endl;

    System sys1;
    sys1.load(path+"2lao.pdb");

    System sys2;
    sys2.load(path+"1viz.pdb");

    Selection sel1(sys1,"name CA");
    cout << "# of CA in 2LAO " << sel1.num() << endl;

    // 1st form of modify
    sel1.modify("resname ALA");
    cout << "# of ALA in 2LAO " << sel1.num() << endl;

    // 2nd form of modify
    sel1.modify(sys2,"name CA");
    cout << "# of CA in 1VIZ " << sel1.num() << endl;

    // Creating empty selection
    Selection sel2;
    cout << "Empty selection: " << sel2.num() << endl;

    cout << "==================================="<<endl;
    cout << "Assignment of selections" << endl;
    cout << "==================================="<<endl;
    sel2 = sel1;
    cout << "After assignment: " << sel2.get_text() << ", " <<sel2.num() << endl;

    cout << "==================================="<<endl;
    cout << " Adding selections to containers"<<endl;
    cout << "==================================="<<endl;
    vector<Selection> sel_list;
    sel_list.push_back(sel1);
    sel_list.push_back(sel2);
    sel_list[1].modify(sys1,"resname THR");

    cout << "vector item 0: " << sel_list[0].get_text() << ", " << sel_list[0].num() << endl;
    cout << "vector item 1: " << sel_list[1].get_text() << ", " << sel_list[1].num() << endl;

    // Making pre-filled container
    cout << "Pre-filled vector:"<<endl;
    vector<Selection> vec2(10);

    cout << "==================================="<<endl;
    cout << " Test residue selection" << endl;
    cout << "==================================="<<endl;
    sel1.each_residue(vec2);
    cout << "Number of residues: " << vec2.size() << endl;
    cout << "Strings and indexes for first 10 residues:" << endl;
    for(int i=0;i<10;++i)
        cout << vec2[i].get_text() << ": " << VectorXi::Map(&vec2[i].get_index()[0], vec2[i].get_index().size()).transpose() << endl;

    cout << "==================================="<<endl;
    cout << " Test IO "<< endl;
    cout << "==================================="<<endl;
    sys1.load(path+"2lao.pdb");
    sys1.load(path+"2lao.pdb");
    Selection all(sys1,"all");
    all.write(path+"test-trj.xtc");

    cout << "Num frames: " << sys1.num_frames() << endl;
    cout << "Num selections: " << sys1.num_selections() << endl;

    // read back
    sys1.load(path+"test-trj.xtc");
    cout << "# frames: " << sys1.num_frames() << endl;
    cout << "# selections: " << sys1.num_selections() << endl;

    Eigen::VectorXf m = sel1.get_mass();
    for(int i=0;i<10;++i) cout << i << " " << m[i] << endl;

    // Test IO of gro files
    System gro(path+"1viz.gro");
    System pdb(path+"1viz.pdb");
    cout << "GRO box:" << endl;
    cout << gro.Box(0) << endl;
    cout << "PDB box:" << endl;
    cout << pdb.Box(0) << endl;
    Selection gro_sel(gro,"all");
    gro_sel.write(path+"res.pdb");
    gro_sel.write(path+"res.gro");

    Vector3f v,a;
    gro.get_box_vectors_angles(0,v,a);
    cout << v.transpose() << " " << a.transpose() << endl;

    cout << "==================================="<<endl;
    cout << " Testing contacts " << endl;
    cout << "==================================="<<endl;
    {

    //System s1("2lao.pdb");
    //Selection dom1(s1,"resid 1-90 or resid 190-238");
    //Selection dom2(s1,"resid 91-189");

    System s1(path+"1viz.pdb");
    Selection dom1(s1,"all");
    Selection dom2(s1,"all");

    vector<Eigen::Vector2i> contacts;
    Grid_searcher(0.5,dom1,contacts,true);
    // Correct number: 481
    ofstream ff((path+"pteros_contacts.dat").c_str());
    for(int i=0;i<contacts.size();++i){
        ff << contacts[i](0) << " " << contacts[i](1) << endl;
    }
    ff.close();

    //Grid_searcher(0.5,dom1,contacts);

    cout << "Found " << contacts.size() << " contacts" << endl;

    // Test periodic searcher
    Grid_searcher(0.5,dom1,contacts,true,true);
    //Grid_searcher(0.5,dom1,contacts,false,true);
    ff.open((path+"pteros_contacts2.dat").c_str());

    for(int i=0;i<contacts.size();++i){
        ff << contacts[i](0) << " " << contacts[i](1) << endl;
    }
    ff.close();


    cout << "Found " << contacts.size() << " periodic contacts" << endl;
    }
    cout << "==================================="<<endl;
    cout << " Testing fitting and rmsd " << endl;
    cout << "==================================="<<endl;

    System op(path+"2lao.pdb");
    System cl(path+"1lst.pdb");

    sel1.modify(op,"name CA and resid 1-238");
    sel2.modify(cl,"name CA and resid 1-238");

    cout << "Before fitting: " << rmsd(sel1,sel2) << endl;

    fit(sel1,sel2);

    cout << "After fitting: " << rmsd(sel1,sel2) << endl;

    sel1.get_name();

*/

/*
    //System_py s;
    cout << "==================================="<<endl;
    cout << " Testing trajectory processor " << endl;
    cout << "==================================="<<endl;

    Processor pr;

    pr.options["trajectory_group"] << string("tr1") << string("tr2") << string("tr3");
    pr.options["topology_file"] << string("dimer_genion.top");
    pr.options["level1"] << 1 << 3.14;
    pr.options["level1/level2"] << 2 << 3;
    pr.options["level1/level2/level3"] << 4;

    //pr.options["trajectory_group/range"] << string("time_range") << 0 << 234;
    cout << "Tree:" << endl;
    ofstream of("js.dat");
    json_spirit::write_stream(pr.options.to_json(),of,true);
    of.close();

    ifstream ff("js.dat");
    mValue vv;
    json_spirit::read_stream(ff,vv);
    ff.close();
    pr.options.from_json(vv);
    cout << "========================\n";
    json_spirit::write_stream(pr.options.to_json(),cout,true);

    string s(pr.options.to_command_line());
    cout << s << endl;
    pr.options.from_command_line(s);
    json_spirit::write_stream(pr.options.to_json(),cout,true);
    cout << "========================\n";

    stringstream ss;
    json_spirit::write_stream(pr.options.to_json(),ss,false);
    cout << ss.str() << endl;
    string s1 = ss.str();
    boost::replace_all(s1,"\"","\\\"");
    pr.options.from_command_line("--json \"" + s1 +"\"");
    json_spirit::write_stream(pr.options.to_json(),cout,true);

    //System sys(path+"dimer_pdb2gmx.gro");
    //pr.set_system(sys);

    //pr.compute();
*/
/*
    cout << "==================================="<<endl;
    cout << " Testing trmsf calculator " << endl;
    cout << "==================================="<<endl;

    RMSF* engine = new RMSF;
    engine->options["trajectory_group"]
    delete engine;
*/

/*
        System sys(path+"/1/dimer_pdb2gmx.gro");
        sys.load(path+"/1/dimer_noPBC_1.xtc",0,10);
        //Selection fitsel(sys,"resid 39 41 42 43 49 51 52 166 170 173 182 184 185 188 214 215 222 223 224 225");
        Selection fitsel(sys,"resid 1");
        Selection all(sys,"all");

        for(int i=5;i<sys.num_frames(); ++i){
            all.set_frame(i);
            Affine3f t = fitsel.fit_transform(i,0);
            cout << t.linear() << endl;
            cout << t.translation().transpose() << endl;
            all.apply_transform(t);
        }
        //all.fit_trajectory(2);

        all.set_frame(0);
        all.write("ref.pdb");
        all.write("aligned.xtc");
*/

    } catch(const Pteros_error& e){ e.print(); }

}

