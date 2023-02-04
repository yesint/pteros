#include "pteros/extras/membrane/lipid_tail.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/utilities.h"
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


LipidTail::LipidTail(LipidTailDescr *descr):
    descr_ptr(descr)
{
    order.resize(descr_ptr->size()-2);
    dihedrals.resize(descr_ptr->size()-3);
}


/**
 * Formulas are taken from:
 * (1) https://pubs.acs.org/doi/10.1021/acs.jctc.7b00643
 * (2) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3882000/
 */
void LipidTail::compute_order_and_dihedrals(const Selection &whole_lipid_sel,
                                            MatrixXf_const_ref normals,
                                            OrderType order_type)
{

    int N = descr_ptr->c_offsets.size();
    bool per_atom_normals = (normals.cols()==1) ? false : true;

    // Compute order
    //atoms:  0 - 1 - 2 - 3 = 4 - 5 - 6
    //bonds:    0   1   2   3   4   5

    const auto& bo = descr_ptr->bond_orders;
    const auto& off = descr_ptr->c_offsets;

    if(order_type == OrderType::SZ){
        // Compute Sz order
        for(int at=1; at<N-1; ++at){
            // Vector from at+1 to at-1
            auto coord1 = whole_lipid_sel.xyz(descr_ptr->c_offsets[at+1]);
            auto coord2 = whole_lipid_sel.xyz(descr_ptr->c_offsets[at-1]);
            auto const& n = per_atom_normals ? normals.col(first_global_index+at) : normals.col(0);
            float ang = angle_between_vectors(coord1-coord2,n);
            order[at-1] = 1.5*pow(cos(ang),2)-0.5;
        }
    } else {
        // Compute deuterium order
        // We iterate over bonds and treat differently single and double bonds
        for(int i=0; i<N-2; ++i){
            if(bo[i]==1){
                // Single bond between atoms i:i+1
                // If next bond is also single, compute order for atom i+1
                if(bo[i+1]==1){
                    // Compute single order for atom i+1
                    /*
                     * C(i)[1]
                     *  \
                     *   C(i+1)[2]---H1,2
                     *  /
                     * C(i+2)[3]
                     */

                    Vector3f p1 = whole_lipid_sel.xyz(off[i]);
                    Vector3f p2 = whole_lipid_sel.xyz(off[i+1]);
                    Vector3f p3 = whole_lipid_sel.xyz(off[i+2]);

                    Vector3f local_z = (p3-p1).normalized();
                    Vector3f local_x = ((p1-p2).cross(p3-p2)).normalized();
                    Vector3f local_y = local_x.cross(local_z);

                    auto const& n = per_atom_normals ? normals.col(first_global_index+i+1) : normals.col(0);
                    float ang_x = angle_between_vectors(local_x,n);
                    float ang_y = angle_between_vectors(local_y,n);
                    float Sxx = 0.5*(3.0*pow(cos(ang_x),2)-1.0);
                    float Syy = 0.5*(3.0*pow(cos(ang_y),2)-1.0);
                    // Instantaneous order
                    // Index in order array is c_atom-1: (i+1)-1=i
                    order[i] = -(2.0*Sxx+Syy)/3.0;
                    // For single bonds there is no difference between ideal
                    // and corrected versions of Scd
                }
                // If next bond is double we do nothing and wait for next iteration
            } else {
                // Double bond between atoms i:i+1
                // Compute double order for atoms i and i+1
                /*
                 * C(i-1)[1]
                 *  \
                 *   C(i)[2]----H1
                 *   ||
                 *   C(i+1)[3]--H2
                 *  /
                 * C(i+2)[4]
                 *
                 * a1 = 0.5*(pi-ang(1,2,3))
                 * a2 = 0.5*(pi-ang(2,3,4))
                 */

                int c1 = off[i-1];
                int c2 = off[i];
                int c3 = off[i+1];
                int c4 = off[i+2];

                Vector3f p1 = whole_lipid_sel.xyz(c1);
                Vector3f p2 = whole_lipid_sel.xyz(c2);
                Vector3f p3 = whole_lipid_sel.xyz(c3);
                Vector3f p4 = whole_lipid_sel.xyz(c4);

                float a1 = 0.5*(M_PI-whole_lipid_sel.angle(c1,c2,c3));
                float a2 = 0.5*(M_PI-whole_lipid_sel.angle(c2,c3,c4));

                // For atom i
                Vector3f local_z = (p3-p2).normalized();
                Vector3f local_x = ((p1-p2).cross(local_z)).normalized();
                Vector3f local_y = local_x.cross(local_z);
                //float ang_x = angle_between_vectors(local_x,normal);
                const auto& n1 = per_atom_normals ? normals.col(first_global_index+i) : normals.col(0);
                float ang_y = angle_between_vectors(local_y,n1);
                float ang_z = angle_between_vectors(local_z,n1);
                float Szz = 0.5*(3.0*pow(cos(ang_z),2)-1.0);
                float Syy = 0.5*(3.0*pow(cos(ang_y),2)-1.0);
                float Syz = 1.5*cos(ang_y)*cos(ang_z);
                // Index in order array is c_atom-1: i-1
                if(order_type == OrderType::SCD_CORR){
                    order[i-1] = -(pow(cos(a1),2)*Syy
                                 + pow(sin(a1),2)*Szz
                                 - 2.0*cos(a1)*sin(a1)*Syz);
                } else { // SCD
                    order[i-1] = -(Szz/4.0 + 3.0*Syy/4.0 - sqrt(3.0)*Syz/2.0);
                }

                // For atom i+1
                //local_z = (p4-p3).normalized();
                local_x = ((p3-p4).cross(local_z)).normalized();
                local_y = local_x.cross(local_z);
                //ang_x = angle_between_vectors(local_x,normal);
                const auto& n2 = per_atom_normals ? normals.col(first_global_index+i+1) : normals.col(0);
                ang_y = angle_between_vectors(local_y,n2);
                ang_z = angle_between_vectors(local_z,n2);
                Szz = 0.5*(3.0*pow(cos(ang_z),2)-1.0);
                Syy = 0.5*(3.0*pow(cos(ang_y),2)-1.0);
                Syz = 1.5*cos(ang_y)*cos(ang_z);
                // Index in order array is c_atom-1: (i+1)-1=i
                if(order_type == OrderType::SCD_CORR){
                    order[i] = -(pow(cos(a2),2)*Syy
                               + pow(sin(a2),2)*Szz
                               + 2.0*cos(a2)*sin(a2)*Syz);
                } else { // SCD
                    order[i] = -(Szz/4.0 + 3.0*Syy/4.0 + sqrt(3.0)*Syz/2.0);
                }
            } // if single/double
        } // for bonds
    } // order_type

    // Compute dihedrals
    for(int at=0; at<N-3; ++at){
        dihedrals[at] = whole_lipid_sel.dihedral(descr_ptr->c_offsets[at],
                                                 descr_ptr->c_offsets[at+1],
                                                 descr_ptr->c_offsets[at+2],
                                                 descr_ptr->c_offsets[at+3],
                                                 noPBC);
    }
}
