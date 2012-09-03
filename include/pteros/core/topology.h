#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <Eigen/Core>
#include <vector>

namespace pteros {

/// Combination rules for non-bond LJ interactions
enum Combination_rule {COMB_RULE_GROMOS, COMB_RULE_AMBER, COMB_RULE_OPLS};

class System;

class Topology
{
    friend class System;
    friend class Selection;
public:
    Topology();
    /// Topology reader
    void read_gromacs_top(std::string fname);
    /// Clear topology
    void clear();
private:
    /// Parent system
    System* system;

    /// Combination rule
    Combination_rule comb_rule;
    /// 1-4 scaling factors
    float scale14_coulomb, scale14_lj;

    /// Bonded parameters
    /// The list of constraints (constrained bonds)
    std::vector<Eigen::Vector2i> bonds;
    /// The list of angles
    std::vector<Eigen::Vector3i> angles;
    /// The list of digedrals
    /// Uses Eigen allocator to ensure correct alignment
    std::vector<Eigen::Vector4i,Eigen::aligned_allocator<Eigen::Vector4i> > dihedrals;
    /// Lenard-Jones params for atoms
    Eigen::VectorXf sigma, epsilon;
    /// Charge
    Eigen::VectorXf charge;
    /// Atom type code
    Eigen::VectorXi type;

    /// Bond distance map.
    /// bond_dist[i,j] shows how many bonds separate atoms i and j.
    Eigen::MatrixXi bond_dist;
};

}

#endif // TOPOLOGY_H
