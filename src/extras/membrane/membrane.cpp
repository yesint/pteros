/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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

#include "pteros/extras/membrane.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "pteros/core/utilities.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <ctime>
#include <fstream>

using namespace std;
using namespace pteros;
using namespace Eigen;



Membrane::Membrane(System *sys, const std::vector<Lipid_descr> &species): system(sys), lipid_species(species)
{
    log = create_logger("membrane");

    Selection all_mid_sel(*system);

    // Creating selections
    for(auto& sp: lipid_species){
        vector<Selection> res;
        system->select(sp.whole_sel_str).split_by_residue(res);
        for(auto& lip: res){
            auto mol = Lipid(lip,sp);
            mol.set_markers();
            all_mid_sel.append(mol.get_mid_marker());
            lipids.push_back(mol);
            resindex_map[mol.get_mid_marker().Resindex(0)] = lipids.size()-1;
        }
    }

    // Compute connectivity    
    std::vector<Selection> leafs;
    all_mid_sel.split_by_connectivity(2.0,leafs,true);
    leaflets.resize(leafs.size());
    leaflets_sel.resize(leafs.size());
    // Set lipid leaflets
    for(int i=0;i<leafs.size();++i){
        leaflets_sel[i].set_system(*system);
        for(auto& a: leafs[i]){
            lipids[resindex_map[a.Resindex()]].set_leaflet(i);
            leaflets[i].push_back(resindex_map[a.Resindex()]);
            leaflets_sel[i].append(a.Index());
        }
    }

    // Unset markers for all lipids
    for(auto& l: lipids) l.unset_markers();

    // Print statictics
    log->info("Number of lipids: {}",lipids.size());
    log->info("Number of leaflets: {}",leafs.size());
}

void Membrane::compute_properties(float d, Vector3f_const_ref external_normal)
{
    bool use_external_normal = false;
    if(external_normal != Vector3f::Zero()) use_external_normal = true;

    // Set markers for all lipids
    for(auto& l: lipids) l.set_markers();

    // Compute per leaflet
    for(int l=0;l<leaflets.size();++l){
        // Get connectivity in this leaflet
        vector<Vector2i> bon;
        search_contacts(d,leaflets_sel[l],bon,false,true);

        // Convert the list of bonds to convenient form
        // atom ==> 1 2 3...
        vector<vector<int> > conn(leaflets_sel[l].size());
        for(int i=0;i<bon.size();++i){
            conn[bon[i](0)].push_back(bon[i](1));
            conn[bon[i](1)].push_back(bon[i](0));
        }

        // Find local normal for each lipid in this leaflet
        for(int i=0;i<leaflets_sel[l].size();++i){

            // Create selection for locality of this lipid
            Selection local = leaflets_sel[l](conn[i]);
            local.append(leaflets_sel[l].Index(i)); // Add central atom
            // Get inertial axes
            Vector3f moments;
            Matrix3f axes;
            local.inertia(moments,axes,true); // Have to use periodic variant

            // Find lipid index corresponding to this mid atom
            int ind = resindex_map[leaflets_sel[l].Resindex(i)];

            Vector3f normal = (use_external_normal) ? external_normal : axes.col(2);

            // Normal is the 3rd axis (shortest)
            // Need to check direction of the normal
            float ang = angle_between_vectors(normal, lipids[ind].get_head_xyz()-lipids[ind].get_tail_xyz());
            if(ang < M_PI_2){
                lipids[ind].normal = normal;
            } else {
                lipids[ind].normal = -normal;
            }
            lipids[ind].normal.normalized();
            lipids[ind].angle = ang;
        }

        //for(int i=0;i<leaflets[l].size();++i){
        //    ind = leaflets[l][i]; // Lipid index
        // }
    }

    // unset markers
    for(auto& l: lipids) l.unset_markers();
}


string tcl_arrow(Vector3f_const_ref p1, Vector3f_const_ref p2, float r, string color){
    stringstream ss;
    Vector3f p = (p2-p1)*0.8+p1;
    ss << "draw color " << color << endl;

    ss << "draw cylinder \"" << p1.transpose()*10.0 << "\" ";
    ss << "\"" << p.transpose()*10.0 << "\" radius " << r << endl;

    ss << "draw cone \"" << p.transpose()*10.0 << "\" ";
    ss << "\"" << p2.transpose()*10.0 << "\" radius " << r*3.0 << endl;

    return ss.str();
}

void Membrane::write_vmd_arrows(const string &fname)
{
    ofstream f(fname);
    for(auto& lip: lipids){
        //vis_dirs << tcl_arrow(markers.XYZ(i),markers.XYZ(i)+props[i].principal_directions.col(0),0.5,"red");
        //vis_dirs << tcl_arrow(markers.XYZ(i),markers.XYZ(i)+props[i].principal_directions.col(1),0.5,"blue");
        f << tcl_arrow(lip.get_mid_xyz(),lip.get_mid_xyz()+lip.normal,0.2,"green");
        f << tcl_arrow(lip.get_tail_xyz(),lip.get_head_xyz(),0.2,"red");
    }
    f.close();

}

Lipid::Lipid(const Selection &sel, const Lipid_descr &descr){
    name = descr.name;
    whole = sel;
    head_sel = whole(descr.head_sel_str);
    tail_sel = whole(descr.tail_sel_str);
    mid_sel = whole(descr.mid_sel_str);
}

void Lipid::set_markers()
{
    // Unwrap this lipid with leading index of position[0]
    whole.unwrap(mid_sel.Index(0)-whole.Index(0));

    // save coords of first atoms
    saved_head0 = head_sel.XYZ(0);
    saved_tail0 = tail_sel.XYZ(0);
    saved_mid0 = mid_sel.XYZ(0);

    // Set markers to COM
    head_sel.XYZ(0) = head_sel.center(true);
    tail_sel.XYZ(0) = tail_sel.center(true);
    mid_sel.XYZ(0) = mid_sel.center(true);
}

void Lipid::unset_markers()
{
    head_sel.XYZ(0) = saved_head0;
    tail_sel.XYZ(0) = saved_tail0;
    mid_sel.XYZ(0) = saved_mid0;
}

