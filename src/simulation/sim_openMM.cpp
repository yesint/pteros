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
#include "OpenMM.h"
#include <vector>
#include <fstream>

class Simulator {
    public:
        Simulator();
        ~Simulator();
        void read_ff(std::string fname);
        void init();
        void step_one();
        void get_state();

    private:
        int Natoms;

        OpenMM::System* sys;
        OpenMM::Context* context;
        OpenMM::Integrator* integrator;
};

using namespace std;
using namespace OpenMM;

Simulator::Simulator(){
    int i,j;
    double v1,v2;

    static bool hasLoadedPlugins = false;
    if (!hasLoadedPlugins) {
        cout << "Loading OpenMM plugins..." << endl;
        vector<string> loadedPlugins =
            OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());
        for (int i = 0; i < loadedPlugins.size(); i++)
            cout << "\tSuccessfully loaded plugin " << loadedPlugins[i] << endl;
        hasLoadedPlugins = true;
    }


    // Create system, context and integrator
    sys = new OpenMM::System();
};

Simulator::~Simulator(){
    delete context;
    delete integrator;
    delete sys;
}


void Simulator::read_ff(std::string fname){
    int i,N;
    double v1,v2,v3,v4,v5,v6;
    int p1,p2,p3,p4;
    string dum, line, str;
    istringstream ss;
    OpenMM::Vec3 vec1,vec2,vec3;

    cout << "Reading FF parameters from " << fname << endl;

    ifstream f(fname.c_str());
    // CMM motion remover
    getline(f,line); ss.clear(); ss.str(line);
    ss >> dum >> p1;

    cout << "\tCM removed once per " << p1 << " steps" << endl;
    if(p1>0) sys->addForce(new CMMotionRemover(p1));

    // Bonds
    getline(f,line); ss.clear(); ss.str(line);
    ss >> dum >> N;
    if(N){
        HarmonicBondForce* bondForce = new HarmonicBondForce();
        sys->addForce(bondForce);
        cout << "\tAdding " << N << " bonds to the system..."<<endl;
        for(i=0;i<N;++i) {
            getline(f,line); ss.clear(); ss.str(line);
            ss >> p1 >> p2 >> v1 >> v2;
            bondForce->addBond(p1,p2,v1,v2);
        }
    }

    // Angles
    getline(f,line); ss.clear(); ss.str(line);
    ss >> dum >> N;
    if(N){
        HarmonicAngleForce* angleForce = new HarmonicAngleForce();
        sys->addForce(angleForce);
        cout << "\tAdding " << N << " angles to the system..."<<endl;
        for(i=0;i<N;++i) {
            getline(f,line); ss.clear(); ss.str(line);
            ss >> p1 >> p2 >> p3 >> v1 >> v2;
            angleForce->addAngle(p1,p2,p3,v1,v2);
        }
    }

    //Proper dihedrals
    getline(f,line); ss.clear(); ss.str(line);
    ss >> dum >> N;
    if(N){
        PeriodicTorsionForce* periodicForce = new PeriodicTorsionForce();
        sys->addForce(periodicForce);
        cout << "\tAdding " << N << " torsion angles to the system..."<<endl;
        for(i=0;i<N;++i) {
            getline(f,line); ss.clear(); ss.str(line);
            ss >> p1 >> p2 >> p3 >> p4 >> v1 >> v2 >> v3;
            periodicForce->addTorsion(p1,p2,p3,p4,v1,v2,v3);
        }
    }

    //RB dihedrals
    getline(f,line); ss.clear(); ss.str(line);
    ss >> dum >> N;
    if(N){
        RBTorsionForce* rbForce = new RBTorsionForce();
        sys->addForce(rbForce);
        cout << "\tAdding " << N << " RB torsion angles to the system..."<<endl;
        for(i=0;i<N;++i) {
            getline(f,line); ss.clear(); ss.str(line);
            ss >> p1 >> p2 >> p3 >> p4 >> v1 >> v2 >> v3 >> v4 >> v5 >> v6;
            rbForce->addTorsion(p1, p2, p3, p4, v1, v2, v3, v4, v5, v6);
        }
    }

    // Non-bonded interactions
    NonbondedForce* nonbondedForce = new NonbondedForce();
    sys->addForce(nonbondedForce);

    // Cut-off type
    getline(f,line); ss.clear(); ss.str(line);
    ss >> dum >> str;
    if(str=="NonbondedForce::CutoffPeriodic"){
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> v1;
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> vec1[0] >> vec1[1] >> vec1[2];
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> vec2[0] >> vec2[1] >> vec2[2];
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> vec3[0] >> vec3[1] >> vec3[2];

        nonbondedForce->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
        nonbondedForce->setCutoffDistance(v1);
        sys->setPeriodicBoxVectors(vec1, vec2, vec3);

        cout << "\tElectrostics type: Periodic cut-off " << v1 << endl;
    } else if(str=="NonbondedForce::CutoffNonPeriodic") {
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> v1;

        nonbondedForce->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
        nonbondedForce->setCutoffDistance(v1);

        cout << "\tElectrostics type: Non-periodic cut-off" << v1 << endl;
    } else if(str=="NonbondedForce::Ewald") {
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> v1;
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> vec1[0] >> vec1[1] >> vec1[2];
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> vec2[0] >> vec2[1] >> vec2[2];
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> vec3[0] >> vec3[1] >> vec3[2];

        nonbondedForce->setNonbondedMethod(NonbondedForce::Ewald);
        nonbondedForce->setCutoffDistance(v1);
        sys->setPeriodicBoxVectors(vec1, vec2, vec3);

        cout << "\tElectrostics type: Ewald; Cut-off = " << v1 << endl;
    } else if(str=="NonbondedForce::PME") {
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> v1;
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> vec1[0] >> vec1[1] >> vec1[2];
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> vec2[0] >> vec2[1] >> vec2[2];
        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> vec3[0] >> vec3[1] >> vec3[2];

        nonbondedForce->setNonbondedMethod(NonbondedForce::PME);
        nonbondedForce->setCutoffDistance(v1);
        sys->setPeriodicBoxVectors(vec1, vec2, vec3);

        cout << "\tElectrostics type: PME; Cut-off = " << v1 << endl;
    }

    //Non-bonds
    getline(f,line); ss.clear(); ss.str(line);
    ss >> dum >> Natoms;
    cout << "\tAdding LJ parameters for " << Natoms << " atoms..." << endl;
    for(i=0;i<Natoms;++i) {
        getline(f,line); ss.clear(); ss.str(line);
        ss >> v1 >> v2 >> v3 >> v4;
        sys->addParticle(v1);
        nonbondedForce->addParticle(v2, v3, v4);
    }

    // Exceptions
    getline(f,line); ss.clear(); ss.str(line);
    ss >> dum >> N;
    cout << "\tAdding " << N << " LJ exceptions..." << endl;
    for(i=0;i<N;++i) {
        getline(f,line); ss.clear(); ss.str(line);
        ss >> p1 >> p2 >> v1 >> v2 >> v3;
        nonbondedForce->addException(p1, p2, v1, v2, v3);
    }

    //GBSA
    getline(f,line); ss.clear(); ss.str(line);
    ss >> dum >> p1;
    if(p1==0){
        cout << "\tGBSA is off" << endl;
    } else {
        GBSAOBCForce* gbsa = new GBSAOBCForce();
        sys->addForce(gbsa);

        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> v1;
        gbsa->setSoluteDielectric(v1);
        cout << "\tGBSA Solute Dielectric: " << v1 << endl;

        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> v2;
        gbsa->setSolventDielectric(v2);
        cout << "\tGBSA Solvent Dielectric: " << v2 << endl;

        getline(f,line); ss.clear(); ss.str(line);
        ss >> dum >> v3;
        cout << "\tGBSA cut-off: " << v3 << endl;
        gbsa->setCutoffDistance(v3);

        cout << "\tAdding GBSA parameters for: " << Natoms << " atoms..." << endl;
        for(i=0;i<Natoms;++i) {
            getline(f,line); ss.clear(); ss.str(line);
            ss >> v1 >> v2 >> v3;
            gbsa->addParticle(v1,v2,v3);
        }
    }

    f.close();
    cout << "Done reading parameters."<<endl;
}

void Simulator::init(){
    cout << "Initializing simulation..." << endl;
    integrator = new VerletIntegrator(0.002);
    context = new Context(*sys, *integrator);

    cout << "Platform: " << context->getPlatform().getName() << endl;

    // Create pteros system
    cout << "\tReading structure file..."<<endl;
    pteros::System pteros_sys("gromacs_files/prot_amber.pdb");
    pteros::Selection all(pteros_sys,"all");

    // Sanity check
    if(all.size()!=Natoms){
        throw pteros::Pteros_error("Inconsistent number of atoms!");
    }

    // Set initial coordinates and velocities
    cout << "\tSetting coordinates and velocities..."<<endl;
    vector<Vec3> pos(Natoms);
    vector<Vec3> vel(Natoms);
    for (int i = 0; i < Natoms; ++i) {
        pos[i] = Vec3(0.1*all.X(i), 0.1*all.Y(i), 0.1*all.Z(i));
        vel[i] = Vec3(0.0, 0.0, 0.0);
    }
    context->setPositions(pos);
    context->setVelocities(vel);

    cout << "Ready to run!" << endl;
}

void Simulator::step_one(){
    integrator->step(1);
}

void Simulator::get_state(){
    int types = 0;
    //types += State::Positions;
    //types += State::Forces;
    types += State::Energy;

    State cur = context->getState(types);

    cout << cur.getPotentialEnergy() << " " << cur.getKineticEnergy() << endl;
}

int main(int argc, char** argv)
{
    clock_t t1,t2;

    Simulator sim;
    sim.read_ff("gromacs_files/OpenMM.params");
    sim.init();

    for(int i=0;i<10;++i){
    //t1 = clock();
    sim.step_one();
    //t2 = clock();
    //cout << "One step takes: " << double(t2-t1)/(double)CLOCKS_PER_SEC << " s." << endl;
    sim.get_state();
    }
}
