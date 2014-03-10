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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <queue>
#include <boost/algorithm/string.hpp> // String algorithms
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include "pteros/core/atom.h"
#include "pteros/core/selection.h"
#include "pteros/core/system.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/grid_search.h"
#include "pteros/core/mol_file.h"
#include "pteros/core/format_recognition.h"
#ifdef USE_POWERSASA
#include "power_sasa.h"
#endif

#include <Eigen/Geometry>
#include <Eigen/Dense>

#include "selection_macro.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

void Selection::allocate_parser(){
    // Parse selection here
    // Parser is heavy object, so if selection is not persistent
    // we will delete it after parsing is complete
    parser.reset(new Selection_parser);
    parser->create_ast(sel_text);    
    parser->apply(system, frame, index);
    if(!parser->has_coord){
        parser.reset();
    }
}

Selection::Selection(){
    system = NULL;
    parser.reset();
};

// Aux function, which creates selection
void Selection::create_internal(System& sys, string& str){
    // Set selection string
    sel_text = str;
    boost::trim(sel_text);

    // Expand macro-definitions in the string
    for(int i=0;i<Nmacro;++i)
        boost::replace_all(sel_text,macro[2*i],macro[2*i+1]);

    // Add selection to sys and make connection
    system = &sys;

    // By default points to frame 0
    frame = 0;

    allocate_parser();
}

// Aux function, which creates selection
void Selection::create_internal(System& sys, int ind1, int ind2){
    // No selection string
    sel_text = "";

    // Add selection to sys and save self-pointer
    system = &sys;

    // By default points to frame 0
    frame = 0;
    // No parser needed
    parser.reset();

    // Populate selection directly
    index.clear();
    for(int i=ind1; i<=ind2; ++i) index.push_back(i);
}


// Aux function for deleting selection
void Selection::delete_internal(){
    // Clear index
    index.clear();
    // If parser is persistent, delete it
    parser.reset();
    // Delete connection
    connection.disconnect();

    sel_text="";
}

// Main constructor
Selection::Selection(System& sys, string str){
    create_internal(sys, str);
}

// Constructor without immediate parsing
Selection::Selection(System& sys){
    // Add selection to sys and save self-pointer
    system = &sys;   
    sel_text = "";
    frame = 0;
    parser.reset();
}

Selection::Selection(System& sys, int ind1, int ind2){
    create_internal(sys, ind1, ind2);
}


// Destructor
Selection::~Selection(){
    // Delete parser if it is persistent
    parser.reset();
    connection.disconnect();
    // All the rest will be destroyed automatically
}

void Selection::append(Selection &sel){
    if(!system) throw Pteros_error("Can't append to undefined system!");
    if(!sel.system) throw Pteros_error("Can't append undefined system!");
    if(sel.system!=system) throw Pteros_error("Can't append atoms from other system!");

    copy(sel.index.begin(),sel.index.end(),back_inserter(index));

    sort(index.begin(),index.end());
    vector<int>::iterator it = unique(index.begin(), index.end());
    index.resize( it - index.begin() );

    sel_text = "";
    parser.reset();
}

void Selection::append(int ind){
    if(!system) throw Pteros_error("Can't append to undefined system!");
    if(ind<0 || ind>=system->num_atoms()) throw Pteros_error("Appended index is out of bonds!");
    index.push_back(ind);
    sel_text = "";
}

// Free memory used by selection.
// Selection is still present in parent.
void Selection::clear(){
    // Clear index
    index.clear();
    // If parser is present (persistent), delete it
    parser.reset();
    sel_text = "";
}

// Modifies both system and text in selection
void Selection::modify(System& sys, string str){
    // If selecton is assigned already, delete and clear it
    if(system){
        delete_internal();
    }
    // re-create selection
    create_internal(sys, str);
}

void Selection::modify(System& sys, int ind1, int ind2){
    // If selecton is assigned already, delete and clear it
    if(system){
        delete_internal();
    }
    // re-create selection
    create_internal(sys, ind1, ind2);
}


// Modify selection if string is provided
void Selection::modify(string str){
    if(system==NULL) Pteros_error("Selection does not belong to the system!");
    sel_text = str;
    boost::trim(sel_text);
    index.clear();
    allocate_parser();
}

void Selection::modify(System &sys){
    if(system){
        delete_internal();
    }
    // Add selection to sys and save self-pointer
    system = &sys;
    frame = 0;
    parser.reset();
    sel_text = "";
}


void Selection::modify(int ind1, int ind2){
    if(system==NULL) Pteros_error("Selection does not belong to the system!");
    parser.reset();
    // re-create selection
    sel_text = "";
    // By default points to frame 0
    frame = 0;
    // Populate selection directly
    index.clear();
    for(int i=ind1; i<=ind2; ++i) index.push_back(i);       
}

void Selection::modify(std::vector<int> &ind){
    if(system==NULL) Pteros_error("Selection does not belong to the system!");
    parser.reset();
    // By default points to frame 0
    frame = 0;
    sel_text = "";
    // Create text and populate selection
    index.clear();
    for(int i=0; i<ind.size(); ++i){
        index.push_back(ind[i]);        
    }
}

void Selection::modify(std::vector<int>::iterator it1, std::vector<int>::iterator it2){
    if(system==NULL) Pteros_error("Selection does not belong to the system!");
    parser.reset();
    // By default points to frame 0
    frame = 0;
    sel_text = "";
    // Populate selection
    index.clear();
    while(it1!=it2){
        index.push_back(*it1);        
        it1++;
    }
}

// Assignment
Selection& Selection::operator=(Selection sel){
    // Sanity check
    if(sel.system==NULL) Pteros_error("Operator '=' with selection, which does not belong to the system!");

    // Kill all current data
    clear();

    // Copy selection text, index and frame
    sel_text = sel.sel_text;
    index = sel.index;
    frame = sel.frame;

    // Add to new parent
    system = sel.system;
    // Connect this selection to notifier if needed
    if(sel.signals_enabled()){
        connection.disconnect();
        connection = system->notify_signal.connect( boost::bind(&Selection::notify_slot,this,_1,_2,_3) );
    }

    if(sel.parser){
        parser.reset(new Selection_parser);
        *parser = *(sel.parser);
    };

    return *this;
}

// Equality operator
bool Selection::operator==(const Selection &other) const {
    return index == other.index;
}

// Copy constructor
Selection::Selection(const Selection& sel){
    if(sel.system==NULL){
        // Making copy of empty selection
        system = NULL;
        sel_text = "";
        parser.reset();
        // And return
        return;
    }

    // Add new data
    sel_text = sel.sel_text;
    index = sel.index;
    frame = sel.frame;
    // Add to new parent
    system = sel.system;
    // Connect this selection to notifier
    if(sel.signals_enabled()){
        connection.disconnect();
        connection = system->notify_signal.connect( boost::bind(&Selection::notify_slot,this,_1,_2,_3) );
    }

    // If parser in sel is persistent, allocate it
    if(sel.parser){
        parser.reset(new Selection_parser);
        *parser = *(sel.parser);
    }
}

// Update selection (re-parse selection text)
void Selection::update(){
    if(sel_text!="") modify(sel_text);
}

// Re-apply AST tree for coordinate-dependent selections
void Selection::apply(){
    if(parser){
        // If parser is persistent, do quick eval using ast tree
        parser->apply(system, frame, index);
    }
}

void Selection::set_frame(int fr){
    if(fr>=0 && fr < system->traj.size()){
        frame = fr;
        // If parser is persistent, do quick update
        // This will only work for coordinate-dependent selections
        apply();
    } else {
        throw Pteros_error("Invalid frame to set!");
    }
}

/////////////////////////
// Get and set functions
/////////////////////////

std::string Selection::get_text() const {
    if(sel_text.size()>0){
        return sel_text;
    } else {
        // If text is empty return dumb indexes
        stringstream ss;
        ss << "index ";
        for(int i=0;i<index.size();++i) ss << index[i] << " ";
        return ss.str();
    }
}

vector<char> Selection::get_chain() const {
    vector<char> res;
    int i,n;
    n = index.size();
    res.resize(n,0);
    for(i=0; i<n; ++i) res[i] = system->atoms[index[i]].chain;
    return res;
}

void Selection::set_chain(const vector<char>& data){
    int i,n;
    n = index.size();
    // Sanity check
    if(data.size()!=n){
        Pteros_error e;
        e << "Invalid data size "<< data.size()
          << " for selection of size " << n;
        throw e;
    }
    for(i=0; i<n; ++i) system->atoms[index[i]].chain = data[i];
}

void Selection::set_chain(char data){
    int i,n;
    n = index.size();
    for(i=0; i<n; ++i) system->atoms[index[i]].chain = data;
}

vector<char> Selection::get_unique_chain() const {
    vector<char> tmp,res;
    int i,n;
    n = index.size();
    tmp.resize(n,0);
    for(i=0; i<n; ++i) tmp[i] = system->atoms[index[i]].chain;
    unique_copy(tmp.begin(),tmp.end(), back_inserter(res));
    return res;
}


vector<int> Selection::get_resid() const {
    vector<int> res;
    int i,n;
    n = index.size();
    res.resize(n,0);
    for(i=0; i<n; ++i) res[i] = system->atoms[index[i]].resid;
    return res;
}

void Selection::set_resid(const vector<int>& data){
    int i,n;
    n = index.size();
    // Sanity check
    if(data.size()!=n){
        Pteros_error e;
        e << "Invalid data size "<< data.size()
          << " for selection of size " << n;
        throw e;
    }
    for(i=0; i<n; ++i) system->atoms[index[i]].resid = data[i];
}

void Selection::set_resid(int data){
    int i,n;
    n = index.size();
    for(i=0; i<n; ++i) system->atoms[index[i]].resid = data;
}

vector<int> Selection::get_resindex() const {
    vector<int> res;
    int i,n;
    n = index.size();
    res.resize(n,0);
    for(i=0; i<n; ++i) res[i] = system->atoms[index[i]].resindex;
    return res;
}

std::vector<float> Selection::get_mass() const {
    vector<float> res;
    int i,n;
    n = index.size();
    res.resize(n);
    for(i=0; i<n; ++i) res[i] = system->atoms[index[i]].mass;
    return res;
}

vector<int> Selection::get_unique_resid() const{
    vector<int> tmp,res;
    int i,n;
    n = index.size();
    tmp.resize(n,0);
    for(i=0; i<n; ++i) tmp[i] = system->atoms[index[i]].resid;
    unique_copy(tmp.begin(),tmp.end(), back_inserter(res));
    return res;
}

vector<int> Selection::get_unique_resindex() const {
    vector<int> tmp,res;
    int i,n;
    n = index.size();
    tmp.resize(n,0);
    for(i=0; i<n; ++i) tmp[i] = system->atoms[index[i]].resindex;
    unique_copy(tmp.begin(),tmp.end(), back_inserter(res));
    return res;
}

void Selection::set_mass(const std::vector<float> m){
    int i,n;
    n = index.size();
    if(m.size()!=n){
        Pteros_error e;
        e << "Invalid data size "<< m.size()
          << " for selection of size " << n;
        throw e;
    }
    for(i=0; i<n; ++i) Mass(i) = m[i];
}

void Selection::set_mass(float data){
    int i,n;
    n = index.size();
    for(i=0; i<n; ++i) system->atoms[index[i]].mass = data;
}


vector<string> Selection::get_name() const {
    vector<string> res;
    int i,n;
    n = index.size();
    res.resize(n);
    for(i=0; i<n; ++i) res[i] = system->atoms[index[i]].name;
    return res;
}


void Selection::set_name(const vector<string>& data){
    int i,n;
    n = index.size();
    // Sanity check
    if(data.size()!=n){
        Pteros_error e;
        e << "Invalid data size "<< data.size()
          << " for selection of size " << n;
        throw e;
    }
    for(i=0; i<n; ++i) system->atoms[index[i]].name = data[i];
}

void Selection::set_name(string& data){
    int i,n;
    n = index.size();
    for(i=0; i<n; ++i) system->atoms[index[i]].name = data;
}

vector<string> Selection::get_resname() const {
    vector<string> res;
    int i,n;
    n = index.size();
    res.resize(n);
    for(i=0; i<n; ++i) res[i] = system->atoms[index[i]].resname;
    return res;
}


void Selection::set_resname(const vector<string>& data){
    int i,n;
    n = index.size();
    // Sanity check
    if(data.size()!=n){
        Pteros_error e;
        e << "Invalid data size "<< data.size()
          << " for selection of size " << n;
        throw e;
    }
    for(i=0; i<n; ++i) system->atoms[index[i]].resname = data[i];
}

void Selection::set_resname(string& data){
    int i,n;
    n = index.size();
    for(i=0; i<n; ++i) system->atoms[index[i]].resname = data;
}


MatrixXf Selection::get_xyz() const {
    MatrixXf res;
    get_xyz(res);
    return res;
}

void Selection::get_xyz(MatrixXf_ref res) const {
    int i,n;
    n = index.size();
    res.resize(3,n);
    for(i=0; i<n; ++i) res.col(i) = system->traj[frame].coord[index[i]];
}


// Compute average structure
MatrixXf Selection::average_structure(int b, int e) const {
    MatrixXf res;
    int i,n,fr;
    n = index.size();
    res.resize(3,n);
    res.fill(0.0);
    if(e==-1) e = system->num_frames()-1;
    if(e<b || b<0 || e>system->num_frames()-1 || e<0) throw Pteros_error("Invalid frame range!");
    cout << "Computing avreage structure from frames: "<<b<<":"<<e<<endl;
    for(fr=b;fr<=e;++fr){
        for(i=0; i<n; ++i) res.col(i) = system->traj[fr].coord[index[i]];
    }
    res /= (e-b+1);
    return res;
}


void Selection::set_xyz(const MatrixXf_cref &coord){
    int n = index.size();
    // Sanity check
    if(coord.cols()!=n){
        Pteros_error e;
        e << "Invalid data size "<< coord.size()
          << " for selection of size " << n;
        throw e;
    }
    for(int i=0; i<n; ++i) XYZ(i) = coord.col(i);
}

vector<float> Selection::get_beta() const{
    vector<float> res;
    int i,n;
    n = index.size();
    res.resize(n);
    for(i=0; i<n; ++i) res[i] = system->atoms[index[i]].beta;
    return res;
}

void Selection::set_beta(std::vector<float>& data){
    int i,n;
    n = index.size();
    for(i=0; i<n; ++i) system->atoms[index[i]].beta = data[i];
}

void Selection::set_beta(float data){
    int i,n;
    n = index.size();
    for(i=0; i<n; ++i) system->atoms[index[i]].beta = data;
}

/*
// Compute bonds using grid search
vector<Vector2i> Selection::compute_bonds(){
    vector<Vector2i> ret;
    Grid_searcher(0.16,*this,ret);
    return ret;
}
*/

////////////////////////////////////////////
// Transformations and inquery functions
////////////////////////////////////////////

// Center of geometry
Vector3f Selection::center(bool mass_weighted, bool periodic) const {
    Vector3f res;
    int i;
    int n = index.size();

    if(n==0) throw Pteros_error("Can't get the center of empty selection!");

    res.fill(0.0);

    if(periodic==false){
        if(mass_weighted){
            float m = 0.0;
            for(i=0; i<n; ++i){
                res += _XYZ(i)*_Mass(i);
                m += _Mass(i);
            }
            return res/m;
        } else {
            for(i=0; i<n; ++i)
                res += _XYZ(i);

            return res/n;
        }
    } else {
        // Periodic center
        // We will find closest periodic images of all points
        // using first point as a reference
        Vector3f ref_point = _XYZ(0);
        if(mass_weighted){
            float m = 0.0;
            for(i=0; i<n; ++i){
                res += system->Box(frame).get_closest_image(_XYZ(i),ref_point) * _Mass(i);
                m += _Mass(i);
            }
            return res/m;
        } else {
            for(i=0; i<n; ++i)
                res += system->Box(frame).get_closest_image(_XYZ(i),ref_point);

            return res/n;
        }
    }
}

// Plain translation
void Selection::translate(const Vector3f_cref &v){
    int i,n = index.size();
    for(i=0; i<n; ++i) XYZ(i) += v;
}


// Rotation with given rotation matrix around 0
void Selection::rotate(const Matrix3f_cref &m){
    int n = index.size();
    for(int i=0; i<n; ++i){
        XYZ(i) = m * XYZ(i);
    }
}

// Rotation around specified axis relative to cm
void Selection::rotate(int axis, float angle){
    if(axis<0 || axis>2) throw Pteros_error("Invalid rotation axis!");
    int n = index.size();
    Affine3f m;
    switch(axis){
    case 0:
        m = AngleAxisf(angle,Vector3f::UnitX());
        break;
    case 1:
        m = AngleAxisf(angle,Vector3f::UnitY());
        break;
    case 2:
        m = AngleAxisf(angle,Vector3f::UnitZ());
        break;
    }
    Vector3f cm = center();
    // Translate to 0
    translate(-cm);
    for(int i=0; i<n; ++i){
        XYZ(i) = m * XYZ(i);
    }
    // Translate back
    translate(cm);
}

// Sequence of rotations around x,y,z relative to zero
void Selection::rotate(const Vector3f_cref &angles, const Vector3f_cref &pivot){
    int n = index.size();
    Affine3f m( AngleAxisf(angles[0],Vector3f::UnitX()) *
               AngleAxisf(angles[1],Vector3f::UnitY()) *
               AngleAxisf(angles[2],Vector3f::UnitZ()) );

    // Translate to pivot
    translate(-pivot);
    for(int i=0; i<n; ++i){
        XYZ(i) = m * XYZ(i);
    }
    // Translate back
    translate(pivot);
}


// Rotation around specified axis relative to given pivot
void Selection::rotate(int axis, float angle, const Vector3f_cref &pivot){
    if(axis<0 || axis>2) throw Pteros_error("Invalid rotation axis!");
    int n = index.size();

    Affine3f m;
    switch(axis){
    case 0:
        m = AngleAxisf(angle,Vector3f::UnitX());
        break;
    case 1:
        m = AngleAxisf(angle,Vector3f::UnitY());
        break;
    case 2:
        m = AngleAxisf(angle,Vector3f::UnitZ());
        break;
    }

    // Translate to 0
    translate(-pivot);
    for(int i=0; i<n; ++i){
        XYZ(i) = m * XYZ(i);
    }
    // Translate back
    translate(pivot);
}

// Rotation around given vector relative to pivot
void Selection::rotate(const Vector3f_cref &direction, float angle, const Vector3f_cref &pivot){
    int n = index.size();
    float ay,az;

    Affine3f m( AngleAxisf(angle,direction.normalized()) );

    // Translate to pivot
    translate(-pivot);
    for(int i=0; i<n; ++i){
        XYZ(i) = m * XYZ(i);
    }
    // Translate back
    translate(pivot);
}

////////////////////////////////////
// RMSD and fitting functions
////////////////////////////////////

// RMSD between two frames
float Selection::rmsd(int fr1, int fr2){
    int n = index.size();
    float res = 0.0;

    if(fr1<0 || fr1>=system->num_frames() ||
            fr2<0 || fr2>=system->num_frames()){
        Pteros_error e;
        e << "RMSD requested for frames" << fr1 << " and "<<fr2
          << " while the valid range is (" << 0<<":"<<system->num_frames()-1<<")";
        throw e;
    }

    for(int i=0; i<n; ++i)
        res += (XYZ(i,fr1)-XYZ(i,fr2)).squaredNorm();

    return sqrt(res/n);
}

//RMSD between current and other frame
float Selection::rmsd(int fr) {
    if(fr<0 || fr>=system->num_frames()){
        Pteros_error e;
        e << "RMSD requested for frame" << fr
          << " while the valid range is (" << 0<<":"<<system->num_frames()-1<<")";
        throw e;
    }
    return rmsd(frame,fr);
}

// Apply transformation
void Selection::apply_transform(Affine3f &t){
    int n = size();

    Vector3f v;
    for(int i=0; i<n; ++i){        
        XYZ(i) = (t * XYZ(i)).eval();
    }
}

Energy_components Selection::non_bond_energy() const
{
    Energy_components e;
    int i,j,n;
    n = size();
    for(i=0; i<n-1; ++i)
        for(j=i+1; j<n; ++j)
            system->add_non_bond_energy(e,i,j,frame,true);
    return e;
}

namespace pteros {

// Non-bond energy between two selections
Energy_components non_bond_energy(Selection& sel1, Selection& sel2, int fr){
    // Check if both selections are from the same system
    if(sel1.get_system()!=sel2.get_system())
        throw Pteros_error("Can't compute non-bond energy between selections from different systems!");

    int n1 = sel1.size();
    int n2 = sel2.size();
    int i,j;

    Energy_components e;

    for(i=0;i<n1;++i)
        for(j=0;i<n2;++j)
            sel1.get_system()->add_non_bond_energy(e,i,j,fr,true);

    return e;
}

// RMSD between two selections (specified frames)
float rmsd(Selection& sel1, int fr1, Selection& sel2, int fr2){
    int n1 = sel1.index.size();
    int n2 = sel2.index.size();
    float res = 0.0;
    if(n1!=n2){
        Pteros_error e;
        e << "Incompatible selections for RMSD of sizes"
          << n1 << "and" << n2;
        throw e;
    }
    if(fr1<0 || fr1>=sel1.system->num_frames() ||
            fr2<0 || fr2>=sel2.system->num_frames()){
        Pteros_error e;
        e << "RMSD requested for frames" << fr1 << " and "<<fr2
          << " while the valid ranges are \n"
          <<"(0:"<<sel1.system->num_frames()-1<<") and "
          <<"(0:"<<sel2.system->num_frames()-1<<")";
        throw e;
    }

    for(int i=0; i<n1; ++i)
        res += (sel1.XYZ(i,fr1)-sel2.XYZ(i,fr2)).squaredNorm();

    return sqrt(res/n1);
}

// RMSD between two selections (current frames)
float rmsd(Selection& sel1, Selection& sel2){
    return rmsd(sel1,sel1.get_frame(),sel2,sel2.get_frame());
}

// Fitting transformation
Affine3f fit_transform(Selection& sel1, Selection& sel2){
    int n1 = sel1.size();
    int n2 = sel2.size();

    Affine3f rot;
    Vector3f cm1, cm2;

    if(n1!=n2){
        Pteros_error e;
        e << "Incompatible selections for fitting of sizes "
          << n1 << "and" << n2;
        throw e;
    }

    // Bring centers to zero
    cm1 = sel1.center(true);
    cm2 = sel2.center(true);
    sel1.translate(-cm1);
    sel2.translate(-cm2);

    // The code below is hacked from GROMACS 3.3.3
    // Used to compute the rotation matrix
    // It computes the rot matrix for two selections, which
    // are centerd at zero.

    int i,j,N,r,c;
    Matrix<float,6,6> omega,om;
    Matrix3f u,vh,vk;
    float xnr,xpc,max_d;
    int nrot;

    omega.fill(0.0);
    om.fill(0.0);
    N = sel1.size();

    //Calculate the matrix U
    u.fill(0.0);
    for(i=0;i<N;++i) // Over atoms in selection
        u += sel1.XYZ(i)*sel2.XYZ(i).transpose()*sel1.Mass(i);

    //Construct omega
    for(r=0; r<6; r++){
        for(c=0; c<=r; c++){
            if (r>=3 && c<3) {
                omega(r,c)=u(r-3,c);
                omega(c,r)=u(r-3,c);
            } else {
                omega(r,c)=0;
                omega(c,r)=0;
            }
        }
    }

    //Finding eigenvalues of omega
    Eigen::SelfAdjointEigenSolver<Matrix<float,6,6> > solver(omega);
    om = solver.eigenvectors();

    /*  Copy only the first two eigenvectors
        The eigenvectors are already sorted ascending by their eigenvalues!
    */
    for(j=0; j<2; j++){
        for(i=0; i<3; i++) {
            vh(j,i)=sqrt(2.0)*om(i,5-j);
            vk(j,i)=sqrt(2.0)*om(i+3,5-j);
        }
    }

    // Calculate the last eigenvector as the cross-product of the first two.
    // This insures that the conformation is not mirrored and
    // prevents problems with completely flat reference structures.

    vh.row(2) = vh.row(0).cross(vh.row(1)) ;
    vk.row(2) = vk.row(0).cross(vk.row(1)) ;

    /* Determine rotational part */
    for(r=0; r<3; r++)
        for(c=0; c<3; c++)
            rot.linear()(c,r) = vk(0,r)*vh(0,c) + vk(1,r)*vh(1,c) + vk(2,r)*vh(2,c);

    // Bring centers back
    sel1.translate(cm1);
    sel2.translate(cm2);

    //Clear translation part. This is important!
    rot.translation().fill(0.0);
    // Add translation part to transform. Note reverse order of translations! This is important.
    return Translation3f(cm2) * rot * Translation3f(-cm1) ;
}

// Fit two selection directly
void fit(Selection& sel1, Selection& sel2){
    Affine3f t = pteros::fit_transform(sel1,sel2);
    sel1.apply_transform(t);
}

} //namespace pteros

// Fit all frames in trajectory
void Selection::fit_trajectory(int ref_frame, int b, int e){
    if(e==-1) e = system->num_frames()-1;
    // Sanity check
    if(b<0 || b>=system->num_frames() || b>e) throw Pteros_error("Invalid frame range!");
    if(ref_frame<0 || ref_frame>=system->num_frames())
        throw Pteros_error("Reference frame is out of range!");

    // Aux reference selection
    Selection ref = *this;
    ref.set_frame(ref_frame);

    Affine3f t;

    int n = system->num_frames();
    for(int fr=0;fr<n;++fr){
        set_frame(fr);
        t = pteros::fit_transform(*this,ref);
        apply_transform(t);
    }
}



// Fitting transformation between two frames of the same selection
Affine3f Selection::fit_transform(int fr1, int fr2){
    // Save current frame
    int cur_frame = get_frame();
    set_frame(fr1); // this points to fr1
    // Create aux selection
    Selection s2(*this);
    s2.set_frame(fr2); // Points to fr2
    // Call fit_transform
    Affine3f t = pteros::fit_transform(*this,s2);
    // Restore frame
    set_frame(cur_frame);
    return t;
}

void Selection::fit(int fr1, int fr2){
    Affine3f t = fit_transform(fr1,fr2);
    apply_transform(t);
}

void Selection::minmax(Vector3f_ref min, Vector3f_ref max) const {
    int i,n,j;
    Vector3f xyz;
    n = index.size();
    min.fill(1e10);
    max.fill(-1e10);
    for(i=0; i<n; ++i){
        xyz = _XYZ(i);
        for(j=0; j<3; ++j){
            if(xyz(j)<min(j)) min(j) = xyz(j);
            if(xyz(j)>max(j)) max(j) = xyz(j);
        }
    }
}

//###############################################
// IO functions
//###############################################

void Selection::write(string fname, int b, int e) {
    if(b<-1 || b>=get_system()->num_frames()) throw Pteros_error("Invalid first frame for writing!");
    if(e<-1 || e>=get_system()->num_frames()) throw Pteros_error("Invalid last frame for writing!");
    if(e<b) throw Pteros_error("Invalid frame range for writing!");
    // -1 is special
    if(b==-1) b=get_frame();
    if(e==-1) e=get_frame();
    cout << "Asked to write frames "<<b<<":"<<e<< endl;

    boost::shared_ptr<Mol_file> f = io_factory(fname,'w');

    if(!f->get_content_type().structure && e!=b){
        throw Pteros_error("Can't write the range of frames to structure file!");
    }    

    cout << "Writing to file '" << fname << "'..." << endl;
    for(int fr=b;fr<=e;++fr){
        set_frame(fr);
        f->write(*this,f->get_content_type());
    }

}

void Selection::each_residue(std::vector<Selection>& sel) const {
    int c,r;
    Selection tmp(*system);
    stringstream ss;

    int n = index.size();

    set<string> m;
    // Cycle over all atoms in this selection and classify them using resindex
    for(int i=0; i<n; ++i){
        ss.str(""); ss.clear();
        ss << "resindex " << system->atoms[index[i]].resindex;
        m.insert(ss.str());
    }
    // Now cycle over this set and make selections
    BOOST_FOREACH(string s, m){
        //cout << s << endl;
        sel.push_back( Selection(*get_system(),s) );
    }
}


MatrixXf Selection::atom_traj(int ind, int b, int e) const {
    if(e==-1) e = system->num_frames()-1;
    // Sanity check
    if(ind<0 || ind>=index.size()) throw Pteros_error("Selection index is out of range!");
    if(b<0 || b>=system->num_frames() || b>e) throw Pteros_error("Invalid frame range!");

    int Nfr = e-b+1;

    MatrixXf ret(3,Nfr);

    for(int fr=b;fr<=e;++fr) ret.col(fr) = system->traj[fr].coord[index[ind]];

    return ret;
}

void Selection::atoms_dup(Selection *res_sel){
    system->atoms_dup(index,res_sel);
}

void Selection::atoms_delete(){
    system->atoms_delete(index);
}

void Selection::distribute(Vector3i_ref ncopies, Vector3f_ref shift){
    Selection tmp(*system);
    Selection res;
    int b,e;
    // Distribute over X
    b = system->num_atoms(); // First added atom
    e = b;
    for(int i = 1; i<ncopies(0); ++i){
        atoms_dup(&res);
        e += index.size()-1;
        res.translate(Vector3f(shift(0)*i,0,0));
    }
    // Distribute over Y
    tmp.modify(b,e);
    tmp.append(*this); // Add this selection to tmp
    for(int i = 1; i<ncopies(1); ++i){
        tmp.atoms_dup(&res);
        e += tmp.size()-1;
        res.translate(Vector3f(0,shift(1)*i,0));
    }
    // Distribute over Z
    tmp.modify(b,e);
    tmp.append(*this); // Add this selection to tmp
    for(int i = 1; i<ncopies(2); ++i){
        tmp.atoms_dup();
        res.translate(Vector3f(0,0,shift(2)*i));
    }
}

void Selection::notify_slot(System_notification code, int b, int e){    
    // Process notifications from the system
    switch(code){
    case TOPOLOGY_CHANGED:
        // Dramatic change! Need to update selection and re-parse selection text
        update();
        break;
    case SYSTEM_CLEARED:
        // Clear selection
        clear();
        // Reset frame to 0
        frame = 0;
        break;
    case FRAMES_DELETED:
        // If we point to one of deleted frames, reset to zero
        if(frame>=b && frame<=e) frame = 0;
        break;
    case COORDINATES_CHANGED:
        // If we point to one of changed frames re-apply selection
        if(frame>=b && frame<=e) apply();
        break;
    case FRAME_CHANGE_REQUESTED:
        frame = b;
    }
}

void Selection::enable_signals(){
    connection.disconnect();
    connection = system->notify_signal.connect( boost::bind(&Selection::notify_slot,this,_1,_2,_3) );
}

void Selection::disable_signals(){
    connection.disconnect();
}

bool Selection::signals_enabled() const {
    return connection.connected();
}

void Selection::split_by_connectivity(float d, std::vector<Selection> &res) {
    //cout << "Splitting by connectivity..." << endl;
    int i,j,k;
    // Clear result
    res.clear();
    // Find all connectivity pairs for given cut-off
    vector<Vector2i> pairs;
    Grid_searcher(d,*this,pairs,false,true);

    // Form a connectivity structure in the form con[i]->1,2,5...
    vector<vector<int> > con(this->size());
    for(i=0; i<pairs.size(); ++i){
        con[pairs[i](0)].push_back(pairs[i](1));
        con[pairs[i](1)].push_back(pairs[i](0));
    }

    // If no connectivity just exit
    if(con.size()==0) return;

    // Mask of already used atoms in sel
    VectorXi mask(size());
    mask.fill(0);

    // Start with first atom in sel and find all connectivity consequtively    
    int ngr = 0;
    while(true){
        // Find first non-used atom, which will be the root for new group
        i = 0;
        while(i<mask.size() && mask(i)>0){ ++i; }
        if(i==mask.size()) break; // All connectivity groups are already found
        // Start new group
        ++ngr;        
        Selection s(*system);
        res.push_back(s);
        res.back().sel_text = "index";
        // Atoms to search
        queue<int> to_search;
        // Add root atom to search set
        to_search.push(i);
        mask(i) = ngr; // mark it as used
        // Now cycle over search set until it exhausts
        while(!to_search.empty()){
            k = to_search.front();
            to_search.pop();
            res.back().index.push_back(_Index(k)); // add it to current selection
            res.back().sel_text += " "+boost::lexical_cast<string>(_Index(k));
            // See all atoms connected to k
            for(int j=0; j<con[k].size(); ++j){
                // if atom is not used, add it to search
                if(mask(con[k][j])==0){
                    to_search.push(con[k][j]);
                    mask(con[k][j]) = ngr; // mark it as used
                }
            }
        }
        // No more atoms to search, the group is ready
    }
}

void Selection::split_by_residue(std::vector<Selection> &res)
{
    // We split selection into several by resindex
    res.clear();
    // Map of resindexes to indexs in selections
    map<int,vector<int> > m;
    for(int i=0; i<size(); ++i){
        m[Resindex(i)].push_back(Index(i));
    }
    // Create selections
    map<int,vector<int> >::iterator it;
    for(it=m.begin();it!=m.end();it++){
        res.push_back(Selection(*system));
        res.back().modify( it->second.begin(), it->second.end() );
    }
}

void Selection::inertia(Vector3f_ref moments, Matrix3f_ref axes, bool periodic) const{
    int n = size();
    int i;
    // Compute the central tensor of inertia. Place it into axes
    axes.fill(0.0);

    Vector3f c = center(true,periodic);

    Vector3f p,d;

    if(periodic){
        for(i=0;i<n;++i){
            // 0 point was used as an anchor in periodic center calculation,
            // so we have to use it as an anchor here as well!
            p = system->Box(frame).get_closest_image(_XYZ(i),_XYZ(0));

            d = p-c;
            axes(0,0) += _Mass(i)*( d(1)*d(1) + d(2)*d(2) );
            axes(1,1) += _Mass(i)*( d(0)*d(0) + d(2)*d(2) );
            axes(2,2) += _Mass(i)*( d(0)*d(0) + d(1)*d(1) );
            axes(0,1) -= _Mass(i)*d(0)*d(1);
            axes(0,2) -= _Mass(i)*d(0)*d(2);
            axes(1,2) -= _Mass(i)*d(1)*d(2);
        }
    } else {
        for(i=0;i<n;++i){
            d = _XYZ(i)-c;
            axes(0,0) += _Mass(i)*( d(1)*d(1) + d(2)*d(2) );
            axes(1,1) += _Mass(i)*( d(0)*d(0) + d(2)*d(2) );
            axes(2,2) += _Mass(i)*( d(0)*d(0) + d(1)*d(1) );
            axes(0,1) -= _Mass(i)*d(0)*d(1);
            axes(0,2) -= _Mass(i)*d(0)*d(2);
            axes(1,2) -= _Mass(i)*d(1)*d(2);
        }
    }
    axes(1,0) = axes(0,1);
    axes(2,0) = axes(0,2);
    axes(2,1) = axes(1,2);

    // Now diagonalize inertia tensor
    Eigen::SelfAdjointEigenSolver<Matrix3f> solver(axes);
    // put output values
    axes = solver.eigenvectors();
    moments = solver.eigenvalues();
}

float Selection::gyration(bool periodic) const {
    int n = size();
    int i;
    float d, a = 0.0, b = 0.0;
    Vector3f c = center(true,periodic);
    for(i=0;i<n;++i){
        if(periodic){
            d = system->Box(frame).distance(_XYZ(i),c);
        } else {
            d = (_XYZ(i)-c).norm();
        }
        a += _Mass(i)*d*d;
        b += _Mass(i);
    }
    return sqrt(a/b);
}

void Selection::wrap(const Vector3i_cref &dims){
    for(int i=0;i<size();++i){
        system->Box(frame).wrap_point(XYZ(i),dims);
    }
}

void Selection::unwrap(const Vector3i_cref &dims){
    Vector3f c = center(true,true);
    for(int i=0;i<size();++i){
        XYZ(i) = system->Box(frame).get_closest_image(XYZ(i),c,true,dims);
    }
}

void Selection::unwrap_bonds(float d, const Vector3i_cref &dims){
    // Find all connectivity pairs for given cut-off
    vector<Vector2i> pairs;
    Grid_searcher(d,*this,pairs,false,true);

    // Form a connectivity structure in the form con[i]->1,2,5...
    vector<vector<int> > con(size());
    for(int i=0; i<pairs.size(); ++i){
        con[pairs[i](0)].push_back(pairs[i](1));
        con[pairs[i](1)].push_back(pairs[i](0));
    }

    // Do all wrapping before    
    wrap();
    //write("0.pdb");

    // Mask of moved atoms
    VectorXi moved(size());
    moved.fill(0);
    int Nmoved = 0;

    // Queue of centers to consider
    set<int> todo;
    // Add first atom to the queue
    todo.insert(0);

    int cur;
    for(;;){
        while(!todo.empty()){
            // Get center from the queue
            cur = *(todo.begin());
            todo.erase(todo.begin()); // And pop it from the queue
            //cout << "(C)" << cur << " " << todo.size() << endl;
            moved[cur] = 1; // Mark as moved
            // Unwrap all atoms bound to cur to position near cur
            for(int i=0; i<con[cur].size(); ++i){
                // We only move atoms, which were not yet moved
                if(moved(con[cur][i])==0){
                    // We don't do wrapping here (passing false) since this will bring atoms back!
                    // We intentially want atoms to be unwrapped                    
                    //Vector3f v = XYZ(con[cur][i]);
                    XYZ(con[cur][i]) = system->Box(frame).get_closest_image(XYZ(con[cur][i]),XYZ(cur),false,dims);
                    /*
                    if((v-XYZ(con[cur][i])).norm()>1e-7){
                        cout << "(0) " << v.transpose() << endl;
                        cout << "(1) " << XYZ(con[cur][i]).transpose() << endl << endl;
                    }
                    */
                    // Add moved atom to centers queue
                    todo.insert(con[cur][i]);
                    //cout << "   (+) " << con[cur][i] << endl;
                    ++Nmoved;                    
                }
            }
        }

        if(Nmoved<size()){
            // Find next not moved and add to the queue
            int i=0;
            while(moved(i)==1){
                ++i;
                if(i>=size()) throw Pteros_error("Wrong connectivity?");
            }
            todo.insert(i);
        } else {
            break; // Done
        }

    }
}

Eigen::Affine3f Selection::principal_transform(bool is_periodic){
    Affine3f rot;
    Vector3f cm = center(true,is_periodic);
    translate(-cm);

    Matrix3f m;
    // Compute principal axes
    Vector3f moments;
    Matrix3f axes;
    inertia(moments,axes,is_periodic);
    // Normalize axes
    for(int i=0;i<3;++i) axes.col(i).normalize();
    // Now orient
    // Step 1. Rotate around Z to move projection of axes(0) to X
    m = AngleAxisf(std::atan2(axes.col(0)(0),axes.col(0)(1)), Vector3f::UnitZ());
    // Step 2. Rotate to superimpose col(0) with X
    m = m* AngleAxisf(std::asin(axes.col(0)(2)), Vector3f::UnitY());
    // Step 3. Apply obtained transform to col(1). This gives new position of col(1)
    // And then bring col(1) to Y
    m = m* AngleAxisf(std::acos( (m*axes.col(1))(1) ), Vector3f::UnitX());

    rot.linear() = m;

    // Bring center back
    translate(cm);

    //Clear translation part. This is important!
    rot.translation().fill(0.0);
    // Add translation part to transform. Note reverse order of translations! This is important.
    return Translation3f(cm) * rot * Translation3f(-cm) ;
}

void Selection::principal_orient(bool is_periodic){
    Affine3f tr = principal_transform(is_periodic);
    apply_transform(tr);
}

#ifdef USE_POWERSASA
float Selection::sasa(float probe_r)
{
    return sasa(probe_r,NULL,NULL,NULL);
}

float Selection::sasa(float probe_r, float *total_volume, vector<float> *area_per_atom, vector<float> *volume_per_atom)
{
    // First obtain the VDW radii of all atoms in selection and add probe radius to them
    vector<float> radii(size());
    for(int i=0; i<size(); ++i) radii[i] = VDW(i) + probe_r;

    bool do_v = total_volume ? true : false;
    bool do_a_per_atom = area_per_atom ? true : false;
    bool do_v_per_atom = volume_per_atom ? true : false;

    // Call POWERSASA
    POWERSASA::PowerSasa<float,Eigen::Vector3f> ps(system->traj[frame].coord, radii, 1, 0, do_v || do_v_per_atom, 0);
    ps.calc_sasa_all();

    float v,surf;

    if(do_v || do_v_per_atom){
        for(int i = 0; i < size(); ++i){
            v = ps.getVol()[i];
            if(do_v_per_atom) (*volume_per_atom)[i] = v;
            if(do_v) *total_volume += v;
        }
    }

    for(int i = 0; i < size(); ++i){
        v = ps.getSasa()[i];
        if(do_a_per_atom) (*area_per_atom)[i] = v;
        surf += v;
    }

    return surf;
}
#endif
