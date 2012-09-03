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

#include "pteros/pteros_core.h"
#include "pteros/simulation/simulation.h"

using namespace std;
using namespace pteros;

////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]){

    string dir("/home/semen/work/Projects/kornelyuk/Sasha/dimer_md/Dimer_htyrn_60kDa/1/");
    string fname(dir+"dimer_pdb2gmx.gro");
    string top(dir+"no_water.top");

    System sys(fname);
    sys.load(dir+"small.xtc");

    Simulation sim;
    sim.load_topology(sys,top);
    sim.params.nb_cutoff = 2.0;
    sim.params.is_periodic = true;
    sim.setup();
    sim.set_frame(1);

    Selection r1(sys,"resid 1");
    Selection r2(sys,"resid 2");

    for(int i=0;i<r1.size();++i) cout << r1.Charge(i) << " " << r1.Mass(i) << endl;
    cout << "-----\n";
    for(int i=0;i<r2.size();++i) cout << r2.Charge(i) << " " << r2.Mass(i) << endl;


    vector<Vector2i> bon;
    Grid_searcher(1.0,r1,r2,bon,true,true);

    //cout << "Updating nlist..."<<endl;
    //sim.update_nlist();
    Energy_components e;
    e = sim.non_bond_energy(bon);
    cout << "Non-bond energy: " << "LJ-14: " << e.lj_14 << "LJ-sr: " << e.lj_sr
        << "q-14: " << e.q_14 << "q-sr: " << e.q_sr << endl;

    return 0;

    //System s("gromacs_files/confout.gro");
    System s("na-cl/conf.gro");
    Selection all(s,"all");

    Simulation ff;
    //ff.load_topology(s,"gromacs_files/processed.top");
    ff.load_topology(s,"na-cl/processed.top");

    //vector<string> v;
    //v.push_back("resid 2");
    //v.push_back("resid 3 to 238");
    //ff.add_excluded_groups(v);

    ff.params.nb_cutoff = 1.0;
    ff.setup();


    cout << "Updating nlist..."<<endl;
    ff.update_nlist();

    cout << "Non-bond energy: " << ff.non_bond_energy().total << endl;

    MatrixXf fval;

    ff.compute_force();
    ff.get_force(fval);

    for(int i=0;i<2;++i){
        cout << i << " " << fval.col(i).transpose() << endl;
    }


    ofstream f("force.dat");
    all.set_frame(0);
    all.write("my-struct.pdb");
    for(int step=0;step<10000;++step){
        //cout << "step " << step << endl;
        ff.verlet_step();
        s.frame_dup(0);
        f <<step*ff.params.time_step << " " << fval.col(0).transpose()
            << " " << fval.col(1).transpose()<<endl;
    }
    all.write("my-traj.xtc");
    f.close();
/*
    Selection sel1(s,"resid 90-190");
    Selection sel2(s,"resid 1-89 or resid 191-238");
    cout << ff.get_energy(sel1,sel2) << endl;

    cout << "2821: " << ff.sigma(2821) << " " << ff.epsilon(2821) << " " << ff.charge(2821) << endl;
    cout << "3646: " << ff.sigma(3646) << " " << ff.epsilon(3646) << " " << ff.charge(3646) << endl;

    // Test forces
    ff.non_bond_force();
    cout << ff.force.transpose() << endl;
    */
}
