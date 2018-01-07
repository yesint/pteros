/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *  
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <boost/algorithm/string.hpp> // String algorithms
#include "pteros/core/atom.h"
#include "pteros/core/selection.h"
#include "pteros/core/system.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "selection_parser.h"
#include "pteros/core/mol_file.h"
// Periodic table from VMD molfile plugins
#include "periodic_table.h"

#ifdef USE_POWERSASA
#include "power_sasa.h"
#endif

#include "sasa.h" // From MDTraj


#include <Eigen/Geometry>
#include <Eigen/Dense>

#include "selection_macro.h"
#include "pteros/core/logging.h"

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

void Selection::sort_and_remove_duplicates()
{
    if(index.size()){
        sort(index.begin(),index.end());
        vector<int>::iterator it = unique(index.begin(), index.end());
        index.resize( it - index.begin() );
        if(index[0]<0) throw Pteros_error("Negative index {} present in Selection!",index[0]);
    } else {
        if(size()==0) LOG()->warn("Selection '{}' is empty! Any call of its methods (except size()) will crash your program!", sel_text);
    }
}

Selection::Selection(){
    system = nullptr;
    parser.reset();
    sel_text = "";
    frame = 0;
};

void expand_macro(string& str){    
    for(int i=0;i<macro.size()/2;++i){
        boost::replace_all(str,macro[2*i].c_str(),macro[2*i+1].c_str());
    }
}

// Main constructor
Selection::Selection(const System &sys, string str, int fr){
    // Set selection string
    sel_text = str;
    boost::trim(sel_text);

    // Expand macro-definitions in the string
    expand_macro(sel_text);

    // Add selection to sys
    system = const_cast<System*>(&sys);

    // point to frame
    frame = fr;

    // Sanity check
    if(frame<0 || frame>=system->num_frames())
        throw Pteros_error("Can't make selection for non-existent frame {} (# of frames is {})", frame, system->num_frames());

    allocate_parser();

    // Show warning if empty selection is created
    if(size()==0) LOG()->warn("Selection '{}' is empty! Any call of its methods (except size()) will crash your program!", sel_text);
}

// Constructor without immediate parsing
Selection::Selection(const System &sys){
    // Add selection to sys
    system = const_cast<System*>(&sys);
    sel_text = "";
    frame = 0;
    parser.reset();
}

Selection::Selection(const System &sys, int ind1, int ind2){
    // No selection string
    sel_text = "";
    // set system
    system = const_cast<System*>(&sys);
    // By default points to frame 0
    frame = 0;
    // No parser needed
    parser.reset();

    if(ind2<ind1) throw Pteros_error("Wrong order of indexes {} and {} (must be ascending)", ind1,ind2);

    // Populate selection directly
    index.reserve(ind2-ind1+1);
    for(int i=ind1; i<=ind2; ++i) index.push_back(i);

    // Show warning if empty selection is created
    if(size()==0)
        LOG()->warn("Selection {}:{} is empty! Any call of its methods (except size()) will crash your program!", ind1, ind2);
}

Selection::Selection(const System &sys, const std::vector<int> &ind){
    // No selection string
    sel_text = "";
    // set system
    system = const_cast<System*>(&sys);
    // By default points to frame 0
    frame = 0;
    // No parser needed
    parser.reset();

    // populate selection    
    index = ind;
    sort_and_remove_duplicates();
}

Selection::Selection(const System &sys, std::vector<int>::iterator it1, std::vector<int>::iterator it2){
    // No selection string
    sel_text = "";
    // set system
    system = const_cast<System*>(&sys);
    // By default points to frame 0
    frame = 0;
    // No parser needed
    parser.reset();

    copy(it1,it2,back_inserter(index));
    sort_and_remove_duplicates();
}

Selection::Selection(const System& sys,
                     const std::function<void (const System &, int, std::vector<int> &)>& callback,
                     int fr)
{
    sel_text = "";
    parser.reset();
    frame = fr;
    // Sanity check
    if(frame<0 || frame>=sys.num_frames())
        throw Pteros_error("Can't make selection for non-existent frame {} (# of frames is {})", frame, sys.num_frames());

    // set system
    system = const_cast<System*>(&sys);
    // call callback
    callback(*system,frame,index);
}

// Destructor
Selection::~Selection(){    
    // Everything (including parser) will be destroyed automatically
}

void Selection::append(const Selection &sel){
    if(!system) throw Pteros_error("Can't append to selection with undefined system!");
    if(!sel.system) throw Pteros_error("Can't append selection with undefined system!");
    if(sel.system!=system) throw Pteros_error("Can't append atoms from other system!");

    copy(sel.index.begin(),sel.index.end(),back_inserter(index));

    sort_and_remove_duplicates();

    sel_text = "";
    parser.reset();
}

void Selection::append(int ind){
    if(!system) throw Pteros_error("Can't append to selection with undefined system!");
    if(ind<0 || ind>=system->num_atoms()) throw Pteros_error("Appended index is out of range!");

    if(find(index.begin(),index.end(),ind)==index.end()){
        index.push_back(ind);
        sort(index.begin(),index.end());
    }

    sel_text = "";
    parser.reset();
}

void Selection::remove(const Selection &sel)
{
    vector<int> tmp;
    set_difference(index.begin(),index.end(),
                   sel.index_begin(),sel.index_end(),back_inserter(tmp));
    index = tmp;
    sel_text = "";
    parser.reset();
}

void Selection::remove(int ind)
{
    vector<int> tmp;
    remove_copy(index.begin(),index.end(),back_inserter(tmp),ind);
    index = tmp;
    sel_text = "";
    parser.reset();
}

void Selection::invert()
{
    // Not textual
    sel_text = "";
    auto old_ind = index;
    index.clear();
    index.reserve(system->num_atoms()-old_ind.size()+1);
    // Assign index
    // We can speed up negation a bit by filling only the "gaps"
    int n = old_ind.size();
    int i,j;

    if(n!=0){

        for(j=0;j<old_ind[0];++j) index.push_back(j); //Before first
        for(i=1;i<n;++i)
            for(j=old_ind[i-1]+1;j<old_ind[i];++j) index.push_back(j); // between any two
        for(j=old_ind[n-1]+1;j<system->num_atoms();++j) index.push_back(j); // after last

    } else {
        for(j=0;j<system->num_atoms();++j) index.push_back(j); //All
    }
}

void Selection::set_system(const System &sys){
    clear();
    system = const_cast<System*>(&sys);
}

// Free memory used by selection.
// Selection is still present in parent.
void Selection::clear(){
    // Clear index
    index.clear();
    // If parser is present (persistent), delete it
    parser.reset();
    sel_text = "";
    frame = 0;
}


Selection Selection::select(string str)
{
    // We create new selection taking full control of its internals
    // We use empty constructor and than repeat actions of sel_text constructor
    // but passing our parent index as subspace
    Selection sub;

    // Set selection string
    sub.sel_text = str;
    boost::trim(sub.sel_text);

    // Expand macro-definitions in the string
    expand_macro(sub.sel_text);

    // Add selection to sys
    sub.system = system;

    // point to frame
    sub.frame = frame;

    // Sanity check
    if(sub.frame<0 || sub.frame>=sub.system->num_frames())        
        throw Pteros_error("Can't make selection for non-existent frame {} (# of frames is {})", sub.frame, sub.system->num_frames());

    // Manually allocate parser with subset from parent index
    sub.parser.reset(new Selection_parser(&index));
    sub.parser->create_ast(sub.sel_text);
    sub.parser->apply(sub.system, sub.frame, sub.index);
    if(!sub.parser->has_coord){
        sub.parser.reset();
    }

    // Show warning if empty selection is created
    if(sub.size()==0)
        LOG()->warn("Selection '{}' is empty! Any call of its methods (except size()) will crash your program!", sub.sel_text);

    // And finally return sub
    return sub;
}

Selection Selection::operator()(string str)
{
    return select(str);
}

Selection Selection::select(int ind1, int ind2)
{
    // ind1 and ind2 are LOCAL indexes, convert them to global and just
    // use normal range constructor
    ind1 = index[ind1];
    ind2 = index[ind2];
    return Selection(*system,ind1,ind2);
}

Selection Selection::operator()(int ind1, int ind2)
{
    return select(ind1,ind2);
}

Selection Selection::select(const std::vector<int> &ind)
{
    // LOCAL indexes are given, convert them to global and just
    // use normal range constructor
    vector<int> glob_ind(ind.size());
    for(int i=0; i<ind.size(); ++i) glob_ind[i] = index[ind[i]];
    return Selection(*system,glob_ind);
}

Selection Selection::operator()(const std::vector<int> &ind)
{
    return select(ind);
}

// Modify selection with new selection string
void Selection::modify(string str, int fr){
    if(system==nullptr) throw Pteros_error("Selection does not belong to any system!");
    sel_text = str;
    boost::trim(sel_text);
    // Expand macro-definitions in the string
    expand_macro(sel_text);

    index.clear();
    frame = fr;
    // Sanity check
    if(frame<0 || frame>=system->num_frames())        
        throw Pteros_error("Can't make selection for non-existent frame {} (# of frames is {})", frame, system->num_frames());
    allocate_parser();    
}

void Selection::modify(int ind1, int ind2){
    if(system==nullptr) throw Pteros_error("Selection does not belong to any system!");
    // no parser needed
    parser.reset();
    // not textual
    sel_text = "";    
    // Populate selection directly
    if(ind2<ind1) throw Pteros_error("Wrong order of indexes {} and {} (must be ascending)", ind1,ind2);
    index.clear();
    for(int i=ind1; i<=ind2; ++i) index.push_back(i);       
}

void Selection::modify(const std::vector<int> &ind){
    if(system==nullptr) throw Pteros_error("Selection does not belong to any system!");
    // no parser needed
    parser.reset();    
    // not textual
    sel_text = "";
    // populate selection    
    index = ind;
    sort_and_remove_duplicates();
}

void Selection::modify(std::vector<int>::iterator it1, std::vector<int>::iterator it2){
    if(system==nullptr) throw Pteros_error("Selection does not belong to any system!");
    // no parser needed
    parser.reset();
    // not textual
    sel_text = "";
    // Populate selection
    index.clear();
    copy(it1,it2,back_inserter(index));
    sort_and_remove_duplicates();
}

void Selection::modify(const std::function<void (const System &, int, std::vector<int> &)>& callback, int fr)
{
    if(system==nullptr) throw Pteros_error("Selection does not belong to any system!");
    // no parser needed
    parser.reset();
    // not textual
    sel_text = "";
    // clear index
    index.clear();
    frame = fr;
    // Sanity check
    if(frame<0 || frame>=system->num_frames())
        throw Pteros_error("Can't make selection for non-existent frame {} (# of frames is {})", frame, system->num_frames());
    // fill
    callback(*system,frame,index);

    sort_and_remove_duplicates();
}

void Selection::modify(const System &sys, string str, int fr){
    set_system(sys);    
    modify(str,fr);
}

void Selection::modify(const System &sys, int ind1, int ind2){
    set_system(sys);
    modify(ind1,ind2);
}

void Selection::modify(const System &sys, const std::vector<int> &ind){
    set_system(sys);
    modify(ind);
}

void Selection::modify(const System &sys, std::vector<int>::iterator it1, std::vector<int>::iterator it2){
    set_system(sys);
    modify(it1,it2);
}

void Selection::modify(const System &sys, const std::function<void (const System &, int, std::vector<int> &)>& callback, int fr)
{
    set_system(sys);
    modify(callback,fr);
}

// Assignment
Selection& Selection::operator=(const Selection& other){
    // Sanity check
    if(other.system==nullptr) Pteros_error("Operator '=' with selection, which does not belong to the system!");

    // Self-assigmnet check
    if(&other==this) return *this;

    // Kill all current data
    clear();

    // Copy selection text, index and frame
    sel_text = other.sel_text;
    index = other.index;
    frame = other.frame;

    // Add to new parent
    system = other.system;

    if(other.parser){
        parser.reset(new Selection_parser);
        *parser = *(other.parser);
    };

    return *this;
}

// Equality operator
bool Selection::operator==(const Selection &other) const {
    return (system == other.system) && (index == other.index);
}

Atom_proxy Selection::operator[](int ind) {
    return Atom_proxy(this,ind);
}

Selection Selection::operator~() const {
    Selection res(*system);
    res.frame = frame;
    // Not textual
    res.sel_text = "";
    // Assign index
    // We can speed up negation a bit by filling only the "gaps"
    int n = index.size();
    int i,j;
    if(n!=0){
        res.index.reserve(system->num_atoms()-n);
        for(j=0;j<index[0];++j) res.index.push_back(j); //Before first
        for(i=1;i<n;++i)
            for(j=index[i-1]+1;j<index[i];++j) res.index.push_back(j); // between any two
        for(j=index[n-1]+1;j<system->num_atoms();++j) res.index.push_back(j); // after last
    } else {
        res.index.resize(system->num_atoms());
        for(j=0;j<system->num_atoms();++j) res.index[j] = j; //All
    }
    return res;
}

namespace pteros {

Selection operator|(const Selection &sel1, const Selection &sel2){
    if(sel1.system!=sel2.system) throw Pteros_error("Can't take logical OR of selections belonging to different systems!");
    if(sel1.frame!=sel2.frame) throw Pteros_error("Can't take logical OR of selections pointing to different frames!");
    // Create resulting selection
    Selection res(*sel1.system);
    // Not text based
    res.sel_text = "";
    // No parser needed
    res.parser.reset();
    // Set frame
    res.frame = sel1.frame;
    // Combine indexes
    std::set_union(sel1.index.begin(),sel1.index.end(),
                   sel2.index.begin(),sel2.index.end(),
                   back_inserter(res.index));
    return res;
}

Selection operator&(const Selection &sel1, const Selection &sel2){
    if(sel1.system!=sel2.system) throw Pteros_error("Can't take logical AND of selections belonging to different systems!");
    if(sel1.frame!=sel2.frame) throw Pteros_error("Can't take logical AND of selections pointing to different frames!");
    // Create resulting selection
    Selection res(*sel1.system);
    // Not text based
    res.sel_text = "";
    // No parser needed
    res.parser.reset();
    // Set frame
    res.frame = sel1.frame;
    // Combine indexes
    std::set_intersection(sel1.index.begin(),sel1.index.end(),
                          sel2.index.begin(),sel2.index.end(),
                          back_inserter(res.index));
    return res;
}


ostream& operator<<(ostream &os, const Selection &sel){
    if(sel.size()>0){
        for(int i=0;i<sel.size()-1;++i) os << sel.Index(i) << " ";
        os << sel.Index(sel.size()-1);
    }
    return os;
}

Selection operator-(const Selection &sel1, const Selection &sel2)
{
    Selection res(*sel1.system);
    res.sel_text = "";
    res.parser.reset();
    res.frame = sel1.frame;
    std::set_difference(sel1.index.begin(),sel1.index.end(),
                        sel2.index.begin(),sel2.index.end(),
                        back_inserter(res.index));
    return res;
}

} // namespace pteros

// Copy constructor
Selection::Selection(const Selection& other){
    if(&other==this) return;

    if(other.system==nullptr){
        // Making copy of empty selection
        system = nullptr;
        sel_text = "";
        parser.reset();
        // And return
        return;
    }

    // Add new data
    sel_text = other.sel_text;
    index = other.index;
    frame = other.frame;
    // Add to new parent
    system = other.system;

    // If parser in sel is persistent, allocate it
    if(other.parser){
        parser.reset(new Selection_parser);
        *parser = *(other.parser);
    }
}

// Update selection (re-parse selection text if it exists)
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
    if(fr<0 || fr >= system->num_frames())
        throw Pteros_error("Invalid frame {} to set! Valid range is 0:", fr, system->num_frames());

    if(frame!=fr){
        frame = fr;
        // If parser is persistent, do quick update
        // This will only work for coordinate-dependent selections
        apply();
    }
}


Selection::iterator Selection::begin(){    
    return iterator(this,0);
}

Selection::iterator Selection::end(){
    return iterator(this,size());
}

/////////////////////////
// Get and set functions
/////////////////////////

#define GET_ATOM_PROP(T,prop) \
vector<T> Selection::get_##prop() const { \
    vector<T> tmp; \
    int i,n; \
    n = index.size(); \
    tmp.resize(n); \
    for(i=0; i<n; ++i) tmp[i] = system->atoms[index[i]].prop; \
    return tmp; \
} \


#define GET_ATOM_PROP_WITH_UNIQUE(T,prop) \
vector<T> Selection::get_##prop(bool unique) const { \
    vector<T> tmp; \
    int i,n; \
    n = index.size(); \
    tmp.resize(n); \
    for(i=0; i<n; ++i) tmp[i] = system->atoms[index[i]].prop; \
    if(unique){ \
        vector<T> res; \
        unique_copy(tmp.begin(),tmp.end(), back_inserter(res)); \
        return res; \
    } \
    return tmp; \
} \


#define SET_ATOM_PROP(T,prop) \
void Selection::set_##prop(const vector<T>& data){ \
    int i,n; \
    n = index.size(); \
    if(data.size()!=n) throw Pteros_error("Invalid data size {} for selection of size {}", data.size(),n); \
    for(i=0; i<n; ++i) system->atoms[index[i]].prop = data[i]; \
} \
void Selection::set_##prop(T data){ \
    int i,n; \
    n = index.size(); \
    for(i=0; i<n; ++i) system->atoms[index[i]].prop = data; \
}


GET_ATOM_PROP_WITH_UNIQUE(char,chain)
SET_ATOM_PROP(char,chain)

GET_ATOM_PROP_WITH_UNIQUE(int,resid)
SET_ATOM_PROP(int,resid)

GET_ATOM_PROP_WITH_UNIQUE(string,name)
SET_ATOM_PROP(string,name)

GET_ATOM_PROP_WITH_UNIQUE(string,resname)
SET_ATOM_PROP(string,resname)

GET_ATOM_PROP_WITH_UNIQUE(int,resindex)

GET_ATOM_PROP(float,mass)
SET_ATOM_PROP(float,mass)

GET_ATOM_PROP(float,beta)
SET_ATOM_PROP(float,beta)

GET_ATOM_PROP(float,occupancy)
SET_ATOM_PROP(float,occupancy)

GET_ATOM_PROP_WITH_UNIQUE(string,tag)
SET_ATOM_PROP(string,tag)


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

MatrixXf Selection::get_xyz(bool make_row_major_matrix) const {
    int n = index.size();
    MatrixXf res;
    if(make_row_major_matrix){
        res.resize(n,3);
        for(int i=0; i<n; ++i) res.row(i) = system->traj[frame].coord[index[i]];
    } else {
        // Column major, default
        res.resize(3,n);
        for(int i=0; i<n; ++i) res.col(i) = system->traj[frame].coord[index[i]];
    }
    return res;
}

void Selection::get_xyz(MatrixXf_ref res) const {
    int i,n;
    n = index.size();
    res.resize(3,n);
    for(i=0; i<n; ++i) res.col(i) = system->traj[frame].coord[index[i]];
}

void Selection::set_xyz(pteros::MatrixXf_const_ref coord){
    int n = index.size();
    // Sanity check
    if(coord.cols()!=n && coord.rows()!=n) throw Pteros_error("Invalid data size {} for selection of size {}", coord.size(),n);
    if(coord.cols()==n){ // Column major, default
        for(int i=0; i<n; ++i) XYZ(i) = coord.col(i);
    } else { // row-major, from python bindings
        for(int i=0; i<n; ++i) XYZ(i) = coord.row(i);
    }
}

// Compute average structure
MatrixXf Selection::average_structure(int b, int e, bool make_row_major_matrix) const {
    MatrixXf res;
    int i,n,fr;
    n = index.size();

    if(e==-1) e = system->num_frames()-1;
    if(e<b || b<0 || e>system->num_frames()-1 || e<0){
        throw Pteros_error("Invalid frame range for average structure!");
    }

    LOG()->debug("Computing avreage structure from frames {}:{}",b,e);

    if(!make_row_major_matrix){
        res.resize(3,n);
        res.fill(0.0);
        for(fr=b;fr<=e;++fr) for(i=0; i<n; ++i) res.col(i) = system->traj[fr].coord[index[i]];
    } else {
        res.resize(n,3);
        res.fill(0.0);
        for(fr=b;fr<=e;++fr) for(i=0; i<n; ++i) res.row(i) = system->traj[fr].coord[index[i]];
    }
    res /= (e-b+1);
    return res;
}


////////////////////////////////////////////
// Transformations and inquery functions
////////////////////////////////////////////

// Center of geometry
Vector3f Selection::center(bool mass_weighted, Array3i_const_ref pbc, int leading_index) const {
    Vector3f res;
    int i;
    int n = index.size();

    if(n==0) throw Pteros_error("Can't get the center of empty selection!");

    res.fill(0.0);

    if( (pbc==0).all() ){
        if(mass_weighted){
            float mass = 0.0;
            #pragma omp parallel
            {                
                Vector3f r(Vector3f::Zero());
                #pragma omp for nowait reduction(+:mass)
                for(i=0; i<n; ++i){
                    r += XYZ(i)*Mass(i);
                    mass += Mass(i);
                }
                #pragma omp critical
                {
                    res += r;
                }
            }
            if(mass==0) throw Pteros_error("Selection has zero mass! Center of mass failed!");
            return res/mass;
        } else {
            #pragma omp parallel
            {
                Vector3f r(Vector3f::Zero()); // local to omp thread
                #pragma omp for nowait
                for(i=0; i<n; ++i) r += XYZ(i);
                #pragma omp critical
                {
                    res += r;
                }
            }

            return res/n;
        }
    } else {
        // Periodic center
        // We will find closest periodic images of all points
        // using first point as a reference
        Vector3f ref_point = XYZ(leading_index);
        Periodic_box& b = system->Box(frame);
        if(mass_weighted){
            float mass = 0.0;
            #pragma omp parallel
            {                
                Vector3f r(Vector3f::Zero());
                #pragma omp for nowait reduction(+:mass)
                for(i=0; i<n; ++i){
                    r += b.get_closest_image(XYZ(i),ref_point,pbc) * Mass(i);
                    mass += Mass(i);
                }
                #pragma omp critical
                {
                    res += r;                    
                }
            }
            if(mass==0) throw Pteros_error("Selection has zero mass! Center of mass failed!");
            return res/mass;
        } else {
            #pragma omp parallel
            {
                Vector3f r(Vector3f::Zero()); // local to omp thread
                #pragma omp for nowait
                for(i=0; i<n; ++i)
                    r += b.get_closest_image(XYZ(i),ref_point,pbc);
                #pragma omp critical
                {
                    res += r;
                }
            }
            return res/n;
        }
    }
}

// Plain translation
void Selection::translate(Vector3f_const_ref v){
    int i,n = index.size();
    #pragma omp parallel for
    for(i=0; i<n; ++i) XYZ(i) += v;
}


// Rotation with given rotation matrix around 0
void Selection::rotate(Matrix3f_const_ref m){
    int n = index.size();
    #pragma omp parallel for
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
    #pragma omp parallel for
    for(int i=0; i<n; ++i){
        XYZ(i) = m * XYZ(i);
    }
    // Translate back
    translate(cm);
}

// Sequence of rotations around x,y,z relative to pivot
void Selection::rotate(Vector3f_const_ref angles, Vector3f_const_ref pivot){
    int n = index.size();
    Affine3f m( AngleAxisf(angles[0],Vector3f::UnitX()) *
               AngleAxisf(angles[1],Vector3f::UnitY()) *
               AngleAxisf(angles[2],Vector3f::UnitZ()) );

    // Translate to pivot
    translate(-pivot);
    #pragma omp parallel for
    for(int i=0; i<n; ++i){
        XYZ(i) = m * XYZ(i);
    }
    // Translate back
    translate(pivot);
}


// Rotation around specified axis relative to given pivot
void Selection::rotate(int axis, float angle, Vector3f_const_ref pivot){
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
    #pragma omp parallel for
    for(int i=0; i<n; ++i){
        XYZ(i) = m * XYZ(i);
    }
    // Translate back
    translate(pivot);
}

// Rotation around given vector relative to pivot
void Selection::rotate(Vector3f_const_ref direction, float angle, Vector3f_const_ref pivot){
    int n = index.size();

    Affine3f m( AngleAxisf(angle,direction.normalized()) );

    // Translate to pivot
    translate(-pivot);
    #pragma omp parallel for
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
float Selection::rmsd(int fr1, int fr2) const{
    int n = index.size();
    float res = 0.0;

    if(fr1<0 || fr1>=system->num_frames() ||
       fr2<0 || fr2>=system->num_frames())
        throw Pteros_error("RMSD requested for frames {}:{} while the valid range is 0:{}", fr1,fr2,system->num_frames()-1);

    #pragma omp parallel for reduction(+:res)
    for(int i=0; i<n; ++i)
        res += (XYZ(i,fr1)-XYZ(i,fr2)).squaredNorm();

    return sqrt(res/n);
}

//RMSD between current and other frame
float Selection::rmsd(int fr) const {
    if(fr<0 || fr>=system->num_frames())
        throw Pteros_error("RMSD requested for frame {} while the valid range is 0:{}", fr,system->num_frames()-1);

    return rmsd(frame,fr);
}

// Apply transformation
void Selection::apply_transform(const Affine3f &t){
    int n = size();    
    #pragma omp parallel for
    for(int i=0; i<n; ++i){        
        XYZ(i) = t * XYZ(i);
    }
}

namespace pteros {

void copy_coord(const Selection &from, int from_fr, Selection &to, int to_fr)
{
    if(from.size()!=to.size())
        throw Pteros_error("Can't copy coordinates between selections of different size!");

    for(int i=0; i<from.size(); ++i){
        to.XYZ(i,to_fr) = from.XYZ(i,from_fr);
    }
}


Vector2f get_energy_for_list(const vector<Vector2i>& pairs, const vector<float>& dist, const System& sys){
    Force_field& ff = const_cast<System&>(sys).get_force_field();
    Vector2f e_total(0,0);
    #pragma omp parallel
    {
        int at1,at2;
        Vector2f e(0,0);
        #pragma omp for nowait
        for(int i=0;i<pairs.size();++i){
            at1 = pairs[i](0);
            at2 = pairs[i](1);
            e += ff.pair_energy(at1, at2, dist[i],
                                sys.Atom_data(at1).charge, sys.Atom_data(at2).charge,
                                sys.Atom_data(at1).type,   sys.Atom_data(at2).type);
        }

        #pragma omp critical
        {
            e_total += e;
        }
    }
    return e_total;
}


Vector2f non_bond_energy(const Selection& sel1,
                         const Selection& sel2,
                         float cutoff,
                         int fr,
                         Array3i_const_ref pbc)
{
    // Check if both selections are from the same system
    if(sel1.get_system()!=sel2.get_system())
        throw Pteros_error("Can't compute non-bond energy between selections from different systems!");

    if(fr<0) fr = sel1.get_frame();

    // Need to set frame fr for both selection to get correct grid search
    int fr1 = sel1.get_frame();
    int fr2 = sel2.get_frame();

    if(fr1!=fr) const_cast<Selection&>(sel1).set_frame(fr);
    if(fr2!=fr) const_cast<Selection&>(sel2).set_frame(fr);

    float d;
    if(cutoff==0){
        d = std::min(sel1.get_system()->get_force_field().rcoulomb, sel1.get_system()->get_force_field().rvdw);
    } else {
        d = cutoff;
    }

    // Perform grid search
    vector<Vector2i> pairs;
    vector<float> dist;
    search_contacts(d,sel1,sel2,pairs,true, (pbc!=0).any() ,&dist);

    // Restore frames if needed
    if(fr1!=fr) const_cast<Selection&>(sel1).set_frame(fr1);
    if(fr2!=fr) const_cast<Selection&>(sel2).set_frame(fr2);

    // Now get energy using pair list and distances
    return get_energy_for_list(pairs,dist,*sel1.get_system());
}


// RMSD between two selections (specified frames)
float rmsd(const Selection& sel1, int fr1, const Selection& sel2, int fr2){
    int n1 = sel1.index.size();
    int n2 = sel2.index.size();
    float res = 0.0;
    if(n1!=n2) throw Pteros_error("Incompatible selections for RMSD of sizes {} and {}", n1,n2);

    if(fr1<0 || fr1>=sel1.system->num_frames() || fr2<0 || fr2>=sel2.system->num_frames())
        throw Pteros_error("RMSD requested for frames {}:{} while the valid range is {}:{}",
                          fr1,fr2,sel1.system->num_frames()-1,sel2.system->num_frames()-1);


    #pragma omp parallel for reduction(+:res)
    for(int i=0; i<n1; ++i)
        res += (sel1.XYZ(i,fr1)-sel2.XYZ(i,fr2)).squaredNorm();

    return sqrt(res/n1);
}

// RMSD between two selections (current frames)
float rmsd(const Selection& sel1, const Selection& sel2){
    return rmsd(sel1,sel1.get_frame(),sel2,sel2.get_frame());
}

// Fitting transformation
Affine3f fit_transform(const Selection& sel1, const Selection& sel2){
    int n1 = sel1.size();
    int n2 = sel2.size();

    if(n1!=n2) throw Pteros_error("Incompatible selections for fitting of sizes {} and {}", n1, n2);

    Affine3f rot;
    Vector3f cm1, cm2;

    // Bring centers to zero
    cm1 = sel1.center(true);
    cm2 = sel2.center(true);

    const_cast<Selection&>(sel1).translate(-cm1);
    const_cast<Selection&>(sel2).translate(-cm2);

    // The code below is hacked from GROMACS 3.3.3
    // Used to compute the rotation matrix
    // It computes the rot matrix for two selections, which
    // are centerd at zero.

    int i,j,N,r,c;
    Matrix<float,6,6> omega,om;
    Matrix3f u,vh,vk;

    omega.fill(0.0);
    om.fill(0.0);
    N = sel1.size();

    //Calculate the matrix U
    u.fill(0.0);
    #pragma omp parallel
    {
        Matrix3f _u(Matrix3f::Zero());
        #pragma omp for nowait
        for(i=0;i<N;++i) // Over atoms in selection
            _u += sel1.XYZ(i)*sel2.XYZ(i).transpose()*sel1.Mass(i);
        #pragma omp critical
        {
            u += _u;
        }
    }

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
    const_cast<Selection&>(sel1).translate(cm1);
    const_cast<Selection&>(sel2).translate(cm2);

    //Clear translation part. This is important!
    rot.translation().fill(0.0);
    // Add translation part to transform. Note reverse order of translations! This is important.
    return Translation3f(cm2) * rot * Translation3f(-cm1) ;
}

// Fit two selection directly
void fit(Selection& sel1, const Selection& sel2){    
    Affine3f t = pteros::fit_transform(sel1,sel2);
    sel1.apply_transform(t);
}

} //namespace pteros


Vector2f Selection::non_bond_energy(float cutoff, Array3i_const_ref pbc) const
{
    float d;
    if(cutoff==0){
        d = std::min(system->get_force_field().rcoulomb, system->get_force_field().rvdw);
    } else {
        d = cutoff;
    }

    // Perform grid search
    vector<Vector2i> pairs;
    vector<float> dist;    
    search_contacts(d,*this,pairs,true, (pbc!=0).any() ,&dist);

    // Now get energy using pair list and distances
    return get_energy_for_list(pairs,dist,*system);
}

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
Affine3f Selection::fit_transform(int fr1, int fr2) const {    
    // Save current frame
    int cur_frame = get_frame();

    const_cast<Selection*>(this)->set_frame(fr1); // this points to fr1
    // Create aux selection
    Selection s2(*this);    
    s2.set_frame(fr2); // Points to fr2

    // Call fit_transform
    Affine3f t = pteros::fit_transform(*this,s2);

    // Restore frame
    const_cast<Selection*>(this)->set_frame(cur_frame);
    return t;
}

void Selection::fit(int fr1, int fr2){
    Affine3f t = fit_transform(fr1,fr2);
    apply_transform(t);
}

void Selection::minmax(Vector3f_ref min, Vector3f_ref max) const {
    int i,n,j;
    Vector3f* p;
    n = index.size();

    float x_min, y_min, z_min, x_max, y_max, z_max;
    x_min = y_min = z_min = 1e10;
    x_max = y_max = z_max = -1e10;

    #pragma omp parallel for reduction(min:x_min,y_min,z_min) reduction(max:x_max,y_max,z_max)
    for(i=0; i<n; ++i){
        p = XYZ_ptr(i);
        if((*p)(0)<x_min) x_min = (*p)(0);
        if((*p)(0)>x_max) x_max = (*p)(0);
        if((*p)(1)<y_min) y_min = (*p)(1);
        if((*p)(1)>y_max) y_max = (*p)(1);
        if((*p)(2)<z_min) z_min = (*p)(2);
        if((*p)(2)>z_max) z_max = (*p)(2);
    }

    min(0) = x_min;
    min(1) = y_min;
    min(2) = z_min;

    max(0) = x_max;
    max(1) = y_max;
    max(2) = z_max;

}

//###############################################
// IO functions
//###############################################

void Selection::write(string fname, int b, int e) {    
    // -1 has special meaning
    if(b==-1) b=get_frame(); // current frame
    if(e==-1) e=system->num_frames()-1; // last frame

    if(b<-1 || b>=get_system()->num_frames()) throw Pteros_error("Invalid first frame for writing!");
    if(e<-1 || e>=get_system()->num_frames()) throw Pteros_error("Invalid last frame for writing!");
    if(e<b) throw Pteros_error("Invalid frame range for writing!");

    auto f = Mol_file::open(fname,'w');

    if(!(f->get_content_type().traj()) && e!=b){
        throw Pteros_error("Can't write the range of frames to structure file!");
    }    

    for(int fr=b;fr<=e;++fr){        
        set_frame(fr);
        f->write(*this,f->get_content_type());
    }
}

void Selection::write(const std::unique_ptr<Mol_file> &handler, Mol_file_content what, int b, int e)
{
    // -1 has special meaning
    if(b==-1) b=get_frame(); // current frame
    if(e==-1) e=system->num_frames()-1; // last frame

    if(b<-1 || b>=get_system()->num_frames()) throw Pteros_error("Invalid first frame for writing!");
    if(e<-1 || e>=get_system()->num_frames()) throw Pteros_error("Invalid last frame for writing!");
    if(e<b) throw Pteros_error("Invalid frame range for writing!");

    if(!(handler->get_content_type().traj()) && e!=b && what.traj()){
        throw Pteros_error("Can't write the range of frames to this file!");
    }

    // First write all except trajectory (if any)
    if(what.atoms() || what.coord()){
        auto c = what;
        c.traj(false);
        handler->write(*this,c);
    }

    // Now write trajectory if asked
    if(what.traj()){
        for(int fr=b;fr<=e;++fr){
            set_frame(fr);
            handler->write(*this, Mol_file_content().traj(true));
        }
    }
}

void Selection::flatten()
{
    parser.reset();
    sel_text = "";
}

string Selection::to_gromacs_ndx(string name)
{
    stringstream s;
    s << "[ " << name << " ]" << endl;
    int n=0;
    for(int ind: index){
        s << ind+1;
        ++n;
        if(n==15){
            n=0;
            s << endl;
        } else {
            s << " ";
        }
    }
    s << endl;
    return s.str();
}


void Selection::each_residue(std::vector<Selection>& sel) const {            
    sel.clear();

    // Resindexes are contigous, so for each atom we search forward and backward to find
    // all atoms of enclosing residue.

    // Set of residues, which are already searched
    set<int> used;

    int b,e; // Begin and end of current residue
    int ind;
    for(int i=0; i<size(); ++i){
        // Skip used
        ind = Resindex(i);
        if(used.count(ind)) continue;
        // Starting global index
        b = e = Index(i);
        // Go backward
        while( b-1>=0 && system->atoms[b-1].resindex == ind){ --b; };
        // Go forward
        while( e+1<system->atoms.size() && system->atoms[e+1].resindex == ind){ ++e; };
        sel.push_back(Selection(*system));

        sel.back().modify(b,e);
        // Mark as used
        used.insert(ind);
    }
}

float Selection::VDW(int ind) const {
    int el = Element_number(ind);
    if(el==0){
        switch(Name(ind)[0]){
            case 'H': return  0.12;
            case 'C': return  0.17;
            case 'N': return  0.155;
            case 'O': return  0.152;
            case 'S': return  0.18;
            case 'P': return  0.18;
            case 'F': return  0.147;
            default:  return  0.15;
        }
    } else {
        // Use periodic table from VMD plugins
        return (el<nr_pte_entries) ? 0.1*pte_vdw_radius[el] : 0.15;
    }
}

string Selection::Element_name(int ind) const {
    int el = Element_number(ind);
    return (el<nr_pte_entries) ? string(pte_label[el]) : "X";
}


MatrixXf Selection::atom_traj(int ind, int b, int e, bool make_row_major_matrix) const {
    if(e==-1) e = system->num_frames()-1;
    // Sanity check
    if(ind<0 || ind>=index.size()) throw Pteros_error("Selection index is out of range!");
    if(b<0 || b>=system->num_frames() || b>e) throw Pteros_error("Invalid frame range!");

    int Nfr = e-b+1;
    MatrixXf ret;

    if(!make_row_major_matrix){
        ret.resize(3,Nfr);
        for(int fr=b;fr<=e;++fr) ret.col(fr) = system->traj[fr].coord[index[ind]];
    } else {
        ret.resize(Nfr,3);
        for(int fr=b;fr<=e;++fr) ret.row(fr) = system->traj[fr].coord[index[ind]];
    }

    return ret;
}

void Selection::split_by_connectivity(float d, std::vector<Selection> &res, bool periodic) {

    vector<vector<int>> con(size());
    res.clear();

    if(d==0){
        // Use bonds from topology
        if(!system->force_field.ready) throw Pteros_error("Can't split by topology: no topology!");
        if(system->force_field.bonds.size()==0) throw Pteros_error("Can't split by topology: no bonds in topology!");

        get_local_bonds_from_topology(con);
    } else {
        // Find all connectivity pairs for given cut-off
        vector<Vector2i> pairs;
        search_contacts(d,*this,pairs,false,periodic);

        // Form a connectivity structure in the form con[i]->1,2,5...
        for(int i=0; i<pairs.size(); ++i){
            con[pairs[i](0)].push_back(pairs[i](1));
            con[pairs[i](1)].push_back(pairs[i](0));
        }
    }

    // Mask of used atoms
    VectorXi used(size());
    used.fill(0);

    // Queue of centers to consider
    set<int> todo;

    // Add first atom to the queue
    int cur = 0;

    todo.insert(cur);
    used[cur] = 1;
    int Nused = 1;

    // Create first selection
    {
        Selection tmp(*system);
        tmp.set_frame(frame);
        tmp.sel_text = "";
        tmp.index.push_back(index[cur]);
        res.push_back(tmp);
    }

    for(;;){
        while(!todo.empty()){
            // Get center from the queue
            cur = *(todo.begin());
            todo.erase(todo.begin()); // And pop it from the queue

            // Find all atoms connected to this one
            for(int i: con[cur]){
                // We only add atoms, which were not yet used as centers
                if(used(i)==0){
                    // Add this atom to selection                    
                    res.back().index.push_back(index[i]);
                    // Add this atom to centers queue
                    todo.insert(i);
                    // Mark as used
                    used[i] = 1;
                    ++Nused;
                }
            }
        }

        if(Nused<size()){
            // Find next not used and add to the queue
            int i=0;
            while(used(i)==1){
                ++i;
                if(i>=size()) throw Pteros_error("UPS! Wrong connectivity?");
            }

            todo.insert(i);
            used[i] = 1;
            Nused++;
            // Start new selection
            Selection tmp(*system);
            tmp.set_frame(frame);
            tmp.sel_text = "";
            tmp.index.push_back(index[i]);
            res.push_back(tmp);
        } else {
            // Done.
            // Indexes of selections are not sorted sort them
            for(auto& s: res) sort(s.index.begin(),s.index.end());
            // And exit
            break;
        }

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

void Selection::split_by_molecule(std::vector<Selection> &res)
{
    if(!system->force_field.ready) throw Pteros_error("Can't split by molecule: no topology!");

    map<int,vector<int>> m;
    int bmol = 0;
    for(int i=0;i<size();++i){
        // Find molecule for this atom
        for(int j=bmol; j<system->force_field.molecules.size(); ++j){
            if(index[i]>=system->force_field.molecules[j](0) && index[i]<=system->force_field.molecules[j](1)){
                bmol=j;
                m[j].push_back(index[i]);
            }
        }
    }

    // Create selections
    map<int,vector<int> >::iterator it;
    for(it=m.begin();it!=m.end();it++){
        res.push_back(Selection(*system));
        res.back().modify( it->second.begin(), it->second.end() );
    }
}

void Selection::split_by_chain(std::vector<Selection> &chains)
{
    // We split selection into several by chain
    chains.clear();
    // Map of resindexes to indexs in selections
    map<char,vector<int> > m;
    for(int i=0; i<size(); ++i){
        m[Chain(i)].push_back(Index(i));
    }
    // Create selections
    map<char,vector<int> >::iterator it;
    for(it=m.begin();it!=m.end();it++){
        chains.push_back(Selection(*system));
        chains.back().modify( it->second.begin(), it->second.end() );
    }
}

void Selection::split_by_contiguous_index(std::vector<Selection> &parts)
{
    parts.clear();
    // Start first contiguous part
    int b = 0, i = 0;
    while(i<size()){
        while(i+1<size() && index[i+1]==index[i]+1) ++i;
        // Part finished
        parts.push_back(Selection(*system));
        parts.back().modify(index[b],index[i]);
        b = i+1;
        i = b;
    }
}

void Selection::split_by_contiguous_residue(std::vector<Selection> &parts)
{
    parts.clear();
    // Start first contiguous part
    int b = 0, i = 0;
    while(i<size()){
        while(i+1<size() && (Resindex(i+1)==Resindex(i)+1 || Resindex(i+1)==Resindex(i)) ) ++i;
        // Part finished
        parts.push_back(Selection(*system));
        parts.back().modify(index[b],index[i]);
        b = i+1;
        i = b;
    }
}

void Selection::inertia(Vector3f_ref moments, Matrix3f_ref axes, Array3i_const_ref pbc, bool leading_index) const{
    int n = size();
    int i;
    // Compute the central tensor of inertia. Place it into axes

    float axes00=0.0, axes11=0.0, axes22=0.0, axes01=0.0, axes02=0.0, axes12=0.0;

    Vector3f c = center(true,pbc,leading_index);

    if( (pbc!=0).any() ){
        Vector3f anchor = XYZ(leading_index);
        Periodic_box& b = system->Box(frame);
        #pragma omp parallel
        {
            Vector3f p,d;
            float m;            
            #pragma omp for reduction(+:axes00,axes11,axes22,axes01,axes02,axes12)
            for(i=0;i<n;++i){
                // 0 point was used as an anchor in periodic center calculation,
                // so we have to use it as an anchor here as well!
                p = b.get_closest_image(XYZ(i),anchor,pbc);
                d = p-c;
                m = Mass(i);
                axes00 += m*( d(1)*d(1) + d(2)*d(2) );
                axes11 += m*( d(0)*d(0) + d(2)*d(2) );
                axes22 += m*( d(0)*d(0) + d(1)*d(1) );
                axes01 += m*d(0)*d(1);
                axes02 += m*d(0)*d(2);
                axes12 += m*d(1)*d(2);
            }
        }
    } else {
        #pragma omp parallel
        {
            Vector3f d;
            float m;
            #pragma omp for reduction(+:axes00,axes11,axes22,axes01,axes02,axes12)
            for(i=0;i<n;++i){
                d = XYZ(i)-c;
                m = Mass(i);
                axes00 += m*( d(1)*d(1) + d(2)*d(2) );
                axes11 += m*( d(0)*d(0) + d(2)*d(2) );
                axes22 += m*( d(0)*d(0) + d(1)*d(1) );
                axes01 += m*d(0)*d(1);
                axes02 += m*d(0)*d(2);
                axes12 += m*d(1)*d(2);
            }
        }
    }
    axes(0,0) = axes00;
    axes(1,1) = axes11;
    axes(2,2) = axes22;
    axes(0,1) = axes(1,0) = -axes01;
    axes(0,2) = axes(2,0) = -axes02;
    axes(1,2) = axes(2,1) = -axes12;

    // Now diagonalize inertia tensor
    Eigen::SelfAdjointEigenSolver<Matrix3f> solver(axes);
    // put output values
    axes = solver.eigenvectors();
    moments = solver.eigenvalues();
}

float Selection::gyration(Array3i_const_ref pbc, bool leading_index) const {
    int n = size();
    int i;
    float d, a = 0.0, b = 0.0;
    Vector3f c = center(true,pbc,leading_index);

    #pragma omp parallel for private(d) reduction(+:a,b)
    for(i=0;i<n;++i){
        if( (pbc!=0).any() ){
            d = system->Box(frame).distance(XYZ(i),c,pbc);
        } else {
            d = (XYZ(i)-c).norm();
        }        
        a += Mass(i)*d*d;
        b += Mass(i);
    }
    return sqrt(a/b);
}

float Selection::distance(int i, int j, Array3i_const_ref pbc) const
{
    return system->distance(Index(i),Index(j),frame, pbc);
}

float Selection::angle(int i, int j, int k, Array3i_const_ref pbc) const
{
    return system->angle(Index(i),Index(j),Index(k),frame,pbc);
}

float Selection::dihedral(int i, int j, int k, int l, Array3i_const_ref pbc) const
{
    return system->dihedral(Index(i),Index(j),Index(k),Index(l),frame,pbc);
}

void Selection::wrap(Array3i_const_ref pbc){
    for(int i=0;i<size();++i){
        Box().wrap_point(XYZ(i),pbc);
    }
}

void Selection::unwrap(Array3i_const_ref pbc, int leading_index){
    Vector3f c;
    if( (pbc==0).all() && leading_index>=0){
        c = XYZ(leading_index);
    } else {
        c = center(true,pbc,leading_index);
    }
    for(int i=0;i<size();++i){
        XYZ(i) = Box().get_closest_image(XYZ(i),c,pbc);
    }
}

int Selection::unwrap_bonds(float d, Array3i_const_ref pbc, int leading_index){
    int Nparts = 1;

    // a connectivity structure in the form con[i]->1,2,5...
    vector<vector<int> > con(size());

    if(d==0){
        // Use bonds from topology
        if(!system->force_field.ready) throw Pteros_error("Can't unwrap by topology: no topology!");
        if(system->force_field.bonds.size()==0) throw Pteros_error("Can't unwrap by topology: no bonds in topology!");

        get_local_bonds_from_topology(con);
    } else {
        // Find all connectivity pairs for given cut-off
        vector<Vector2i> pairs;
        search_contacts(d,*this,pairs,false,true); // Periodic by definition

        for(int i=0; i<pairs.size(); ++i){
            con[pairs[i](0)].push_back(pairs[i](1));
            con[pairs[i](1)].push_back(pairs[i](0));
        }
    }

    // Mask of used atoms
    VectorXi used(size());
    used.fill(0);

    // Queue of centers to consider
    set<int> todo;

    // Leading atom coordinates
    Vector3f leading;

    todo.insert(leading_index);
    used[leading_index] = 1;
    int Nused = 1;
    Periodic_box& b = system->Box(frame);

    for(;;){
        while(!todo.empty()){
            // Get center from the queue
            int cur = *(todo.begin());
            todo.erase(todo.begin()); // And pop it from the queue
            leading = XYZ(cur);

            // Find all atoms connected to this one
            for(int i=0; i<con[cur].size(); ++i){
                // We only add atoms, which were not yet used as centers
                if(used(con[cur][i])==0){
                    // Unwrap atom
                    XYZ(con[cur][i]) = b.get_closest_image(XYZ(con[cur][i]),leading,pbc);
                    // Add this atom to centers queue
                    todo.insert(con[cur][i]);
                    // Mark as used
                    used[con[cur][i]] = 1;
                    ++Nused;
                }
            }
        }

        if(Nused<size()){
            // Find next not used atom and add to the queue
            int i=0;
            while(used(i)==1){
                ++i;
                if(i>=size()) throw Pteros_error("Wrong connectivity?");
            }

            todo.insert(i);
            used[i] = 1;
            Nused++;

            Nparts++;
        } else {
            // Done. Exit loop.
            break;
        }

    }

    return Nparts;
}

Eigen::Affine3f Selection::principal_transform(Array3i_const_ref pbc, bool leading_index) const {
    Affine3f rot;
    Vector3f cm = center(true,pbc,leading_index);
    // We need to const-cast in order to call non-const translate() from const method
    // It's Ok because we'll translate back later
    const_cast<Selection*>(this)->translate(-cm);

    Matrix3f m;
    // Compute principal axes
    Vector3f moments;
    Matrix3f axes;
    inertia(moments,axes,pbc,leading_index);
    // Normalize axes
    for(int i=0;i<3;++i) axes.col(i).normalize();
    // Now orient
    // Step 1. Rotate around Z to move projection of col(0) to X
    m = AngleAxisf(std::atan2(axes.col(0)(0),axes.col(0)(1)), Vector3f::UnitZ());
    // Step 2. Rotate to superimpose col(0) with X
    m = m* AngleAxisf(std::asin(axes.col(0)(2)), Vector3f::UnitY());
    // Step 3. Apply obtained transform to col(1). This gives new position of col(1)
    // And then bring col(1) to Y
    m = m* AngleAxisf(std::acos( (m*axes.col(1))(1) ), Vector3f::UnitX());

    rot.linear() = m;

    // Bring center back
    // We need to const-cast in order to call non-const translate() from const method
    const_cast<Selection*>(this)->translate(cm);

    //Clear translation part. This is important!
    rot.translation().fill(0.0);
    // Add translation part to transform. Note reverse order of translations! This is important.
    return Translation3f(cm) * rot * Translation3f(-cm) ;
}

void Selection::principal_orient(Array3i_const_ref pbc, bool leading_index){
    Affine3f tr = principal_transform(pbc,leading_index);
    apply_transform(tr);
}

#ifdef USE_POWERSASA

//=====================================================================================
// aux classes which allows to use Selection as container for coordinates with iterator

template<class ValueType>
class Selection_container_it_t {
public:
    typedef ValueType value_type;
    typedef int difference_type;
    typedef ValueType* pointer;
    typedef ValueType& reference;
    typedef std::forward_iterator_tag iterator_category;

    Selection_container_it_t(Selection* sel, int n) {parent = sel; pos = n;}

    Selection_container_it_t operator++(int junk) { Selection_container_it_t tmp = *this; ++pos; return tmp; }
    Selection_container_it_t& operator++() { ++pos; return *this; }
    reference operator*() const { return parent->XYZ(pos); }
    pointer operator->() { return parent->XYZ_ptr(pos); }
    bool operator==(const Selection_container_it_t& rhs) { return pos == rhs.pos && parent == rhs.parent; }
    bool operator!=(const Selection_container_it_t& rhs) { return pos != rhs.pos || parent != rhs.parent; }

    operator Selection_container_it_t<const Eigen::Vector3f>() const { return *this; }
private:
    Selection* parent;
    int pos;
};


class Selection_coord_container {
public:
    Selection_coord_container(Selection& sel): parent(&sel){}

    typedef Selection_container_it_t<Eigen::Vector3f> iterator;
    typedef Selection_container_it_t<const Eigen::Vector3f> const_iterator;

    iterator begin(){ return iterator(parent,0); }
    const_iterator begin() const{ return const_iterator(parent,0); }
    const_iterator cbegin() const{ return const_iterator(parent,0); }
    iterator end(){ return iterator(parent,parent->size()); }
    const_iterator end() const { return const_iterator(parent,parent->size()); }
    const_iterator cend() const { return const_iterator(parent,parent->size()); }

    int size() const {return parent->size();}
private:
    Selection* parent;
};

//==============================================================================


float Selection::powersasa(float probe_r, vector<float> *area_per_atom,
                           float *total_volume, vector<float> *volume_per_atom) const
{
    // First obtain the VDW radii of all atoms in selection and add probe radius to them
    vector<float> radii(size());
    for(int i=0; i<size(); ++i) radii[i] = VDW(i) + probe_r;

    bool do_v = total_volume ? true : false;
    bool do_a_per_atom = area_per_atom ? true : false;
    bool do_v_per_atom = volume_per_atom ? true : false;

    // Create aux container which allows not to copy selection coordinates locally
    Selection_coord_container cont(const_cast<Selection&>(*this));

    // Call POWERSASA
    POWERSASA::PowerSasa<float,Eigen::Vector3f> ps(cont, radii, 1, 0, do_v || do_v_per_atom, 0);
    ps.calc_sasa_all();

    float v,surf;

    if(do_v || do_v_per_atom){
        if(volume_per_atom) volume_per_atom->resize(size());
        for(int i = 0; i < size(); ++i){
            v = ps.getVol()[i];
            if(do_v_per_atom) (*volume_per_atom)[i] = v;
            if(do_v) *total_volume += v;
        }
    }

    if(area_per_atom) area_per_atom->resize(size());
    for(int i = 0; i < size(); ++i){        
        v = ps.getSasa()[i];
        cout << i << " " << v << endl;
        if(do_a_per_atom) (*area_per_atom)[i] = v;
        surf += v;
    }

    return surf;
}

#else
float Selection::powersasa(float probe_r, float *total_volume,
                      vector<float> *area_per_atom, vector<float> *volume_per_atom) const
{
    throw Pteros_error("Pteros is compiled without powersasa support!");
}

#endif


// sasa implementation from MDTraj

float Selection::sasa(float probe_r, vector<float> *area_per_atom, int n_sphere_points) const
{
    // SASA code works for a range of frames. We use only current frame here

    // Make contigous array of coordinates
    vector<Vector3f> coord(size());
    for(int i=0;i<size();++i) coord[i] = XYZ(i);

    vector<float> radii(size());
    for(int i=0; i<size(); ++i) radii[i] = VDW(i) + probe_r;

    vector<int> atom_mapping(size());
    for(int i=0; i<size(); ++i) atom_mapping[i]=i;

    vector<float> out;
    float* out_ptr;

    if(area_per_atom){
        area_per_atom->resize(size());
        out_ptr = area_per_atom->data();
    } else {
        out.resize(size());
        out_ptr = out.data();
    }

    for(int i=0; i<size(); ++i) out_ptr[i]=0.0;

    ::sasa(1,size(),(float*)coord.data(),radii.data(),n_sphere_points,atom_mapping.data(),size(),out_ptr);

    return Map<VectorXf>(out_ptr,size()).sum();
}


void Selection::get_local_bonds_from_topology(vector<vector<int>>& con){
    int bind = Index(0);
    int eind = Index(size()-1);
    int a1,a2;
    auto bit = std::begin(index);
    auto cur_b = bit;
    auto eit = std::end(index);
    vector<int>::iterator it1,it2;
    for(int i=0;i<system->force_field.bonds.size();++i){
        a1 = system->force_field.bonds[i](0);
        a2 = system->force_field.bonds[i](1);
        if(a1>=bind && a1<=eind && a2>=bind && a2<=eind){
            it1 = std::find(cur_b,eit,a1);
            cur_b = it1;
            it2 = std::find(cur_b,eit,a2);
            con[it1-bit].push_back(it2-bit);
            con[it2-bit].push_back(it1-bit);
        }
    }
}

