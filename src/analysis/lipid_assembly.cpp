/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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


#include "pteros/analysis/lipid_assembly.h"
#include "pteros/core/pteros_error.h"
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace pteros;
using namespace Eigen;

Lipid_assembly::Lipid_assembly(Selection &sel, std::string head_marker_atom, float d){
    create(sel,head_marker_atom,d);
}

void Lipid_assembly::create(Selection &sel, std::string head_marker_atom, float d, float bilayer_cutoff){
    source_ptr = &sel;
    System* sys = sel.get_system();
    //Selection all(*sys,"all");
    // Create selection from marker atoms
    Selection markers(*sys,"("+sel.get_text()+") and "+head_marker_atom);    
    int N = markers.size();

    // Find connctivity of marker atoms
    vector<Vector2i> bon;
    Grid_searcher(d,markers,bon,false,true);

    // Convert the list of bonds to convenient form
    vector<vector<int> > conn(N);
    for(int i=0;i<bon.size();++i){
        conn[bon[i](0)].push_back(bon[i](1));
        conn[bon[i](1)].push_back(bon[i](0));
    }

    // Now cycle over connectivity map    
    Vector3f moments;
    Matrix3f axes;
    int atom;
    Selection local(*sys);          

    surface_normals.resize(markers.size());
    surface_mean_angle.resize(markers.size());
    surface_curvature.resize(markers.size());

    for(int i=0;i<markers.size();++i){
        //cout << "(>>) " << i << endl;
        // Selection for all neighbours including itself
        local.clear();
        for(int j=0; j<conn[i].size(); ++j){
            local.append(markers.Index(conn[i][j]));
        }
        local.append(markers.Index(i)); // append atom itself

        // Compute the inertia axes        
        local.inertia(moments,axes,true);        
        // axes.col(2) is a local normal
        surface_normals[i] = axes.col(2);
    }

    // Compute mean angle between the normals of local
    float a,ang;
    float c,curv;
    for(int i=0;i<markers.size();++i){
        //cout << "(**) " << i << endl;
        ang = 0;
        curv = 0;
        for(int j=0; j<conn[i].size(); ++j){
            // Angle between conn[i][j] and i
            a = acos(surface_normals[i].dot(surface_normals[conn[i][j]])/(surface_normals[i].norm() * surface_normals[conn[i][j]].norm()));
            if(a>M_PI_2) a = M_PI-a;
            ang += a;

            c = 2.0*sin(a) / sys->distance(markers.XYZ(i),markers.XYZ(conn[i][j]),0,true);
            curv += c;

        }
        if(conn[i].size()>0){
            ang /= conn[i].size();
            curv /= conn[i].size();
        }

        markers.Occupancy(i) = curv;

        surface_curvature[i] = curv;
        surface_mean_angle[i] = ang;
    }

    // All markers with occupancy < cutoff correspond to bilayer part
    // Obtain monolayers
    Selection bi(*sys,"occupancy < " + boost::lexical_cast<string>(bilayer_cutoff));

    markers.write("colored.pdb");
}

