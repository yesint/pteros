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

void Lipid_assembly::create(Selection &sel, std::string head_marker_atom, float d){
    source_ptr = &sel;
    System* sys = sel.get_system();
    Selection all(*sys,"all");
    // Create selection from marker atoms
    Selection markers(*sys,"("+sel.get_text()+") and "+head_marker_atom);    
    int N = markers.size();
    markers.set_beta(1);
    markers.set_mass(1); // CG atoms may have no mass assigned

    markers.unwrap_bonds(d,Vector3i(1,0,1));
    all.write("wrapped.pdb");

    sys->Box(0).col(0) *= 2.0;
    sys->Box(0).col(2) *= 2.0;

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


    for(int i=0;i<1;++i){
        // Selection for all neighbours including itself
        local.clear();
        for(int j=0; j<conn[i].size(); ++j){
            local.append(markers.Index(conn[i][j]));
        }
        local.append(markers.Index(i)); // append atom itself

        // Compute the inertia axes        
        local.inertia(moments,axes,true);        
        // axes.col(2) is a local normal

        // Select a cylinder along the normal
        local.set_beta(0); // To exclude local set

        cout << "draw cylinder \""
             << (markers.XYZ(i)-20*axes.col(2)).transpose()*10 << "\" \""
             << (markers.XYZ(i)+20*axes.col(2)).transpose()*10
             << "\" radius 1" << endl;

        stringstream ss;
        ss << head_marker_atom
           << " and beta > 0 "
           << " and distance vector "
           << markers.XYZ(i).transpose() << " " << axes.col(2).transpose()
           << " < " << d;
        Selection other(*sys,ss.str());
        cout << other.get_text() << endl;

        cout << "===========" << endl;
        for(int k=0;k<local.size();++k) cout << local.Index(k) << " ";
        cout << endl;
        for(int k=0;k<other.size();++k) cout << other.Index(k) << " ";
        cout << endl;
        cout << "===========" << endl;

        cout << local.center().transpose()*10  << " " << other.center().transpose()*10 << endl;

        // Get the distance between local and other
        float dist = sys->distance(local.center(true,true),other.center(true,true),0);
        cout << i << " " << dist << endl;
    }



    cout << "**********" << endl;
    Selection s(*sys,"name PO4 and distance vector 13.54 3.38 25.76 -1.12 0 1  < 2");
    for(int k=0;k<s.size();++k) cout << s.Index(k) << " ";
    cout << endl;


}

