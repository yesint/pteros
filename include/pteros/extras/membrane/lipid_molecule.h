#pragma once

#include "pteros/core/selection.h"
#include "lipid_tail.h"
//#include "lipid_species.h"
#include "local_patch.h"
#include "quad_surface.h"
//#include "pteros/extras/lipid_membrane.h"

namespace pteros {

class LipidMembrane;
class LipidSpecies;

class LipidMolecule {
    friend class LipidMembrane;
    friend class PerSpeciesProperties;
    friend class LipidGroup;
public:
    LipidMolecule(const Selection& lip_mol, LipidSpecies *sp, int ind, LipidMembrane* parent);

    void add_to_group(int gr);

    Selection whole_sel;
    Selection head_marker_sel;
    Selection tail_marker_sel;
    Selection surf_marker_sel;
    Eigen::Vector3f tail_head_vector;

    // Tails
    std::vector<LipidTail> tails;

    // Instanteneous properties of lipid
    Eigen::Vector3f normal;
    float tilt;
    float gaussian_curvature;
    float mean_curvature;
    float mono_thickness;

    // Properties from Voronoi tesselation
    int coord_number;
    float area;
    std::vector<int> neib;

    // Inclusion neibours
    std::vector<int> inclusion_neib;

    // Fiting resutls
    Eigen::Vector3f smoothed_surf_marker;
    float quad_fit_rms;

private:

    int id;
    LipidMembrane* membr_ptr;
    LipidSpecies* species_ptr;

    // Set markers to current COM coordinates of marker seletions
    void set_markers();    

    // Coordinates of markers
    Eigen::Vector3f head_marker, tail_marker, surf_marker;    

    // Staff related to local patch computations
    LocalPatch patch;
    // Quadric surface
    QuadSurface surf;
};

}
