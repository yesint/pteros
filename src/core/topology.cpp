#include "pteros/core/topology.h"

using namespace pteros;

Topology::Topology()
{
}

void Topology::clear(){
    scale14_coulomb = scale14_lj = 1.0;
    bonds.clear();
    angles.resize(0);
    dihedrals.resize(0);
    sigma.resize(0);
    epsilon.resize(0);
    charge.resize(0);
    type.resize(0);
    bond_dist.resize(0,0);
}
