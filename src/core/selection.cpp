/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
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


#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <regex>
#include "pteros/core/atom.h"
#include "pteros/core/selection.h"
#include "pteros/core/system.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "selection_parser.h"
#include "pteros/core/file_handler.h"
#include "pteros/core/utilities.h"

#ifdef USE_POWERSASA
#include "power_sasa.h"
#endif

#include "sasa.h" // From MDTraj

#include <Eigen/Geometry>
#include <Eigen/Dense>

#include "selection_macro.h"
#include "pteros/core/logging.h"

// DSSP
#include "pteros_dssp_wrapper.h"


using namespace std;
using namespace pteros;
using namespace Eigen;


void Selection::allocate_parser(){
    // Parse selection here
    // Parser is heavy object, so if selection is not persistent
    // we will delete it after parsing is complete
    parser.reset(new SelectionParser);
    parser->create_ast(sel_text,system);
    parser->apply_ast(frame, _index);
    if(!parser->has_coord){
        parser.reset();
    }
}

void Selection::sort_and_remove_duplicates()
{
    if(_index.size()){
        sort(_index.begin(),_index.end());
        vector<int>::iterator it = unique(_index.begin(), _index.end());
        _index.resize( it - _index.begin() );
        if(_index[0]<0) throw PterosError("Negative index {} present in Selection!",_index[0]);
    } else {
        LOG()->debug("Selection '{}' is empty! Any call of its methods (except size()) will crash your program!", sel_text);
    }
}

Selection::Selection(){
    system = nullptr;
    parser.reset();
    sel_text = "";
    frame = 0;
};

void expand_macro(string& str){    
    for(int i=0;i<selection_macro.size()/2;++i){
        std::regex e(selection_macro[2*i]);
        str = std::regex_replace(str,e,selection_macro[2*i+1]);
    }
}

// Main constructor
Selection::Selection(const System &sys, string str, int fr){
    // Set selection string
    sel_text = str;
    str_trim_in_place(sel_text);

    // Expand macro-definitions in the string
    expand_macro(sel_text);

    // Add selection to sys
    system = const_cast<System*>(&sys);

    // point to frame
    frame = fr;

    // Sanity check
    if(frame<0 || frame>=system->num_frames())
        throw PterosError("Can't make selection for non-existent frame {} (# of frames is {})", frame, system->num_frames());

    allocate_parser();

    // Show warning if empty selection is created
    if(size()==0) LOG()->debug("Selection '{}' is empty! Any call of its methods (except size()) will crash your program!", sel_text);
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

    if(ind2<ind1) throw PterosError("Wrong order of indexes {} and {} (must be ascending)", ind1,ind2);

    // Populate selection directly
    _index.reserve(ind2-ind1+1);
    for(int i=ind1; i<=ind2; ++i) _index.push_back(i);

    // Show warning if empty selection is created
    if(size()==0)
        LOG()->debug("Selection {}:{} is empty! Any call of its methods (except size()) will crash your program!", ind1, ind2);
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
    _index = ind;
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

    copy(it1,it2,back_inserter(_index));
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
        throw PterosError("Can't make selection for non-existent frame {} (# of frames is {})", frame, sys.num_frames());

    // set system
    system = const_cast<System*>(&sys);
    // call callback
    callback(*system,frame,_index);
}

// Destructor
Selection::~Selection(){    
    // Everything (including parser) will be destroyed automatically
}

void Selection::append(const Selection &sel){
    if(!system) throw PterosError("Can't append to selection with undefined system!");
    if(!sel.system) throw PterosError("Can't append selection with undefined system!");
    if(sel.system!=system) throw PterosError("Can't append atoms from other system!");

    copy(sel._index.begin(),sel._index.end(),back_inserter(_index));

    sort_and_remove_duplicates();

    sel_text = "";
    parser.reset();
}

void Selection::append(int ind){
    if(!system) throw PterosError("Can't append to selection with undefined system!");
    if(ind<0 || ind>=system->num_atoms()) throw PterosError("Appended index is out of range!");

    if(find(_index.begin(),_index.end(),ind)==_index.end()){
        _index.push_back(ind);
        sort(_index.begin(),_index.end());
    }

    sel_text = "";
    parser.reset();
}

void Selection::append(const std::vector<Selection> &sel_vec)
{
    if(!system) throw PterosError("Can't append to selection with undefined system!");
    for(const auto& sel: sel_vec){
        if(!sel.system) throw PterosError("Can't append selection with undefined system!");
        if(sel.system!=system) throw PterosError("Can't append atoms from other system!");
    }

    for(const auto& sel: sel_vec){
        copy(sel._index.begin(),sel._index.end(),back_inserter(_index));
    }

    sort_and_remove_duplicates();

    sel_text = "";
    parser.reset();
}

void Selection::remove(const Selection &sel)
{
    vector<int> tmp;
    set_difference(_index.begin(),_index.end(),
                   sel.index_begin(),sel.index_end(),back_inserter(tmp));
    _index = tmp;
    sel_text = "";
    parser.reset();
}

void Selection::remove(int ind)
{
    vector<int> tmp;
    remove_copy(_index.begin(),_index.end(),back_inserter(tmp),ind);
    _index = tmp;
    sel_text = "";
    parser.reset();
}

void Selection::invert()
{
    // Not textual
    sel_text = "";
    auto old_ind = _index;
    _index.clear();
    _index.reserve(system->num_atoms()-old_ind.size()+1);
    // Assign index
    // We can speed up negation a bit by filling only the "gaps"
    int n = old_ind.size();
    int i,j;

    if(n!=0){

        for(j=0;j<old_ind[0];++j) _index.push_back(j); //Before first
        for(i=1;i<n;++i)
            for(j=old_ind[i-1]+1;j<old_ind[i];++j) _index.push_back(j); // between any two
        for(j=old_ind[n-1]+1;j<system->num_atoms();++j) _index.push_back(j); // after last

    } else {
        for(j=0;j<system->num_atoms();++j) _index.push_back(j); //All
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
    _index.clear();
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
    str_trim_in_place(sub.sel_text);

    // Expand macro-definitions in the string
    expand_macro(sub.sel_text);

    // Add selection to sys
    sub.system = system;

    // point to frame
    sub.frame = frame;

    // Sanity check
    if(sub.frame<0 || sub.frame>=sub.system->num_frames())        
        throw PterosError("Can't make selection for non-existent frame {} (# of frames is {})", sub.frame, sub.system->num_frames());

    // Manually allocate parser with subset from parent index
    sub.parser.reset(new SelectionParser(&_index));
    sub.parser->create_ast(sub.sel_text,sub.system);
    sub.parser->apply_ast(sub.frame, sub._index);
    if(!sub.parser->has_coord){
        sub.parser.reset();
    }

    // Show warning if empty selection is created
    if(sub.size()==0)
        LOG()->debug("Selection '{}' is empty! Any call of its methods (except size()) will crash your program!", sub.sel_text);

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
    ind1 = _index[ind1];
    ind2 = _index[ind2];
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
    for(int i=0; i<ind.size(); ++i) glob_ind[i] = _index[ind[i]];
    return Selection(*system,glob_ind);
}

Selection Selection::operator()(const std::vector<int> &ind)
{
    return select(ind);
}

// Modify selection with new selection string
void Selection::modify(string str, int fr){
    if(system==nullptr) throw PterosError("Selection does not belong to any system!");
    sel_text = str;
    str_trim_in_place(sel_text);
    // Expand macro-definitions in the string
    expand_macro(sel_text);

    _index.clear();
    frame = fr;
    // Sanity check
    if(frame<0 || frame>=system->num_frames())        
        throw PterosError("Can't make selection for non-existent frame {} (# of frames is {})", frame, system->num_frames());
    allocate_parser();    
}

void Selection::modify(int ind1, int ind2){
    if(system==nullptr) throw PterosError("Selection does not belong to any system!");
    // no parser needed
    parser.reset();
    // not textual
    sel_text = "";    
    // Populate selection directly
    if(ind2<ind1) throw PterosError("Wrong order of indexes {} and {} (must be ascending)", ind1,ind2);
    _index.clear();
    _index.reserve(ind2-ind1+1);
    for(int i=ind1; i<=ind2; ++i) _index.push_back(i);
}

void Selection::modify(const std::vector<int> &ind){
    if(system==nullptr) throw PterosError("Selection does not belong to any system!");
    // no parser needed
    parser.reset();    
    // not textual
    sel_text = "";
    // populate selection    
    _index = ind;
    sort_and_remove_duplicates();
}

void Selection::modify(std::vector<int>::iterator it1, std::vector<int>::iterator it2){
    if(system==nullptr) throw PterosError("Selection does not belong to any system!");
    // no parser needed
    parser.reset();
    // not textual
    sel_text = "";
    // Populate selection
    _index.clear();
    copy(it1,it2,back_inserter(_index));
    sort_and_remove_duplicates();
}

void Selection::modify(const std::function<void (const System &, int, std::vector<int> &)>& callback, int fr)
{
    if(system==nullptr) throw PterosError("Selection does not belong to any system!");
    // no parser needed
    parser.reset();
    // not textual
    sel_text = "";
    // clear index
    _index.clear();
    frame = fr;
    // Sanity check
    if(frame<0 || frame>=system->num_frames())
        throw PterosError("Can't make selection for non-existent frame {} (# of frames is {})", frame, system->num_frames());
    // fill
    callback(*system,frame,_index);

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
    if(other.system==nullptr) PterosError("Operator '=' with selection, which does not belong to the system!");

    // Self-assigmnet check
    if(&other==this) return *this;

    // Kill all current data
    clear();

    // Copy selection text, index and frame
    sel_text = other.sel_text;
    _index = other._index;
    frame = other.frame;

    // Add to new parent
    system = other.system;

    if(other.parser){
        parser.reset(new SelectionParser);
        *parser = *(other.parser);
    };

    return *this;
}

// Equality operator
bool Selection::operator==(const Selection &other) const {
    return (system == other.system) && (_index == other._index);
}

AtomHandler Selection::operator[](int ind) {
    return AtomHandler(*system,_index[ind],frame);
}

AtomHandler Selection::operator[](const std::pair<int, int> &ind_fr) {
    return AtomHandler(*system,_index[ind_fr.first],ind_fr.second);
}

Selection Selection::operator~() const {
    Selection res(*system);
    res.frame = frame;
    // Not textual
    res.sel_text = "";
    // Assign index
    // We can speed up negation a bit by filling only the "gaps"
    int n = _index.size();
    int i,j;
    if(n!=0){
        res._index.reserve(system->num_atoms()-n);
        for(j=0;j<_index[0];++j) res._index.push_back(j); //Before first
        for(i=1;i<n;++i)
            for(j=_index[i-1]+1;j<_index[i];++j) res._index.push_back(j); // between any two
        for(j=_index[n-1]+1;j<system->num_atoms();++j) res._index.push_back(j); // after last
    } else {
        res._index.resize(system->num_atoms());
        for(j=0;j<system->num_atoms();++j) res._index[j] = j; //All
    }
    return res;
}

namespace pteros {

Selection operator|(const Selection &sel1, const Selection &sel2){
    if(sel1.system!=sel2.system) throw PterosError("Can't take logical OR of selections belonging to different systems!");
    if(sel1.frame!=sel2.frame) throw PterosError("Can't take logical OR of selections pointing to different frames!");
    // Create resulting selection
    Selection res(*sel1.system);
    // Not text based
    res.sel_text = "";
    // No parser needed
    res.parser.reset();
    // Set frame
    res.frame = sel1.frame;
    // Combine indexes
    std::set_union(sel1._index.begin(),sel1._index.end(),
                   sel2._index.begin(),sel2._index.end(),
                   back_inserter(res._index));
    return res;
}

Selection operator&(const Selection &sel1, const Selection &sel2){
    if(sel1.system!=sel2.system) throw PterosError("Can't take logical AND of selections belonging to different systems!");
    if(sel1.frame!=sel2.frame) throw PterosError("Can't take logical AND of selections pointing to different frames!");
    // Create resulting selection
    Selection res(*sel1.system);
    // Not text based
    res.sel_text = "";
    // No parser needed
    res.parser.reset();
    // Set frame
    res.frame = sel1.frame;
    // Combine indexes
    std::set_intersection(sel1._index.begin(),sel1._index.end(),
                          sel2._index.begin(),sel2._index.end(),
                          back_inserter(res._index));
    return res;
}


ostream& operator<<(ostream &os, const Selection &sel){
    if(sel.size()>0){
        for(int i=0;i<sel.size()-1;++i) os << sel.index(i) << " ";
        os << sel.index(sel.size()-1);
    }
    return os;
}

Selection operator-(const Selection &sel1, const Selection &sel2)
{
    Selection res(*sel1.system);
    res.sel_text = "";
    res.parser.reset();
    res.frame = sel1.frame;
    std::set_difference(sel1._index.begin(),sel1._index.end(),
                        sel2._index.begin(),sel2._index.end(),
                        back_inserter(res._index));
    return res;
}

bool check_selection_overlap(const std::vector<Selection> &sel_vec)
{
    vector<int> inters;
    for (int i=0; i<sel_vec.size()-1; ++i){
        for (int j=i+1; j<sel_vec.size(); ++j){
            if(sel_vec[i].size()>0 && sel_vec[j].size()>0){
                set_intersection(sel_vec[i].index_begin(), sel_vec[i].index_end(),
                                 sel_vec[j].index_begin(), sel_vec[j].index_end(),
                                 back_inserter(inters));
            }
            if (!inters.empty()) return true;
        }
    }
    return false;
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
    _index = other._index;
    frame = other.frame;
    // Add to new parent
    system = other.system;

    // If parser in sel is persistent, allocate it
    if(other.parser){
        //parser.reset(new Selection_parser);
        //*parser = *(other.parser);
        update(); // Needed to update internal references in parse AST
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
        parser->apply_ast(frame, _index);
    }
}

void Selection::set_frame(int fr){
    if(fr<0 || fr >= system->num_frames())
        throw PterosError("Invalid frame {} to set! Valid range is 0:", fr, system->num_frames());

    if(frame!=fr){
        frame = fr;
        // If parser is persistent, do quick update
        // This will only work for coordinate-dependent selections
        apply();
    }
}


Selection::iterator Selection::begin(){    
    return iterator(*this,0);
}

Selection::iterator Selection::end(){
    return iterator(*this,size());
}

/////////////////////////
// Get and set functions
/////////////////////////

#define GET_ATOM_PROP(T,prop) \
vector<T> Selection::get_##prop() const { \
    vector<T> tmp; \
    int i,n; \
    n = _index.size(); \
    tmp.resize(n); \
    for(i=0; i<n; ++i) tmp[i] = system->atoms[_index[i]].prop; \
    return tmp; \
} \


#define GET_ATOM_PROP_WITH_UNIQUE(T,prop) \
vector<T> Selection::get_##prop(bool unique) const { \
    vector<T> tmp; \
    int i,n; \
    n = _index.size(); \
    tmp.resize(n); \
    for(i=0; i<n; ++i) tmp[i] = system->atoms[_index[i]].prop; \
    if(unique){ \
        vector<T> res; \
        sort(tmp.begin(),tmp.end()); \
        unique_copy(tmp.begin(),tmp.end(), back_inserter(res)); \
        return res; \
    } \
    return tmp; \
} \


#define SET_ATOM_PROP(T,prop) \
void Selection::set_##prop(const vector<T>& data){ \
    int i,n; \
    n = _index.size(); \
    if(data.size()!=n) throw PterosError("Invalid data size {} for selection of size {}", data.size(),n); \
    for(i=0; i<n; ++i) system->atoms[_index[i]].prop = data[i]; \
} \
void Selection::set_##prop(T data){ \
    int i,n; \
    n = _index.size(); \
    for(i=0; i<n; ++i) system->atoms[_index[i]].prop = data; \
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

GET_ATOM_PROP(float,charge)
SET_ATOM_PROP(float,charge)


GET_ATOM_PROP_WITH_UNIQUE(string,tag)
SET_ATOM_PROP(string,tag)


std::string Selection::get_text() const {
    if(sel_text.size()>0){
        return sel_text;
    } else {
        // If text is empty return dumb indexes
        stringstream ss;
        ss << "index ";
        for(int i=0;i<_index.size();++i) ss << _index[i] << " ";
        return ss.str();
    }
}

MatrixXf Selection::get_xyz(bool make_row_major_matrix) const {
    int n = _index.size();
    MatrixXf res;
    if(make_row_major_matrix){
        res.resize(n,3);
        for(int i=0; i<n; ++i) res.row(i) = system->traj[frame].coord[_index[i]];
    } else {
        // Column major, default
        res.resize(3,n);
        for(int i=0; i<n; ++i) res.col(i) = system->traj[frame].coord[_index[i]];
    }
    return res;
}

void Selection::set_xyz(pteros::MatrixXf_const_ref coord){
    int n = _index.size();
    // Sanity check
    if(coord.cols()!=n && coord.rows()!=n) throw PterosError("Invalid data size {} for selection of size {}", coord.size(),n);
    if(coord.cols()==n){ // Column major, default
        for(int i=0; i<n; ++i) xyz(i) = coord.col(i);
    } else { // row-major, from python bindings
        for(int i=0; i<n; ++i) xyz(i) = coord.row(i);
    }
}

MatrixXf Selection::get_vel(bool make_row_major_matrix) const {
    if(!system->traj[frame].has_vel()) throw PterosError("System has no velocities");

    int n = _index.size();
    MatrixXf res;
    if(make_row_major_matrix){
        res.resize(n,3);
        for(int i=0; i<n; ++i) res.row(i) = system->traj[frame].vel[_index[i]];
    } else {
        // Column major, default
        res.resize(3,n);
        for(int i=0; i<n; ++i) res.col(i) = system->traj[frame].vel[_index[i]];
    }
    return res;
}

void Selection::set_vel(pteros::MatrixXf_const_ref data){
    int n = _index.size();
    // Sanity check
    if(data.cols()!=n && data.rows()!=n) throw PterosError("Invalid data size {} for selection of size {}", data.size(),n);

    system->traj[frame].vel.resize(n);

    if(data.cols()==n){ // Column major, default
        for(int i=0; i<n; ++i) vel(i) = data.col(i);
    } else { // row-major, from python bindings
        for(int i=0; i<n; ++i) vel(i) = data.row(i);
    }
}


MatrixXf Selection::get_force(bool make_row_major_matrix) const {
    if(!system->traj[frame].has_force()) throw PterosError("System has no forces");

    int n = _index.size();
    MatrixXf res;
    if(make_row_major_matrix){
        res.resize(n,3);
        for(int i=0; i<n; ++i) res.row(i) = system->traj[frame].force[_index[i]];
    } else {
        // Column major, default
        res.resize(3,n);
        for(int i=0; i<n; ++i) res.col(i) = system->traj[frame].force[_index[i]];
    }
    return res;
}

void Selection::set_force(pteros::MatrixXf_const_ref data){
    int n = _index.size();
    // Sanity check
    if(data.cols()!=n && data.rows()!=n) throw PterosError("Invalid data size {} for selection of size {}", data.size(),n);

    system->traj[frame].force.resize(n);

    if(data.cols()==n){ // Column major, default
        for(int i=0; i<n; ++i) force(i) = data.col(i);
    } else { // row-major, from python bindings
        for(int i=0; i<n; ++i) force(i) = data.row(i);
    }
}


float Selection::get_total_charge() const {
    float q = 0.0;
    for(int i=0; i<size(); ++i) q += system->atoms[_index[i]].charge;
    return q;
}

// Compute average structure
MatrixXf Selection::average_structure(int b, int e, bool make_row_major_matrix) const {
    MatrixXf res;
    int i,n,fr;
    n = _index.size();

    if(e==-1) e = system->num_frames()-1;
    if(e<b || b<0 || e>system->num_frames()-1 || e<0){
        throw PterosError("Invalid frame range for average structure!");
    }

    LOG()->debug("Computing avreage structure from frames {}:{}",b,e);

    if(!make_row_major_matrix){
        res.resize(3,n);
        res.fill(0.0);
        for(fr=b;fr<=e;++fr) for(i=0; i<n; ++i) res.col(i) = system->traj[fr].coord[_index[i]];
    } else {
        res.resize(n,3);
        res.fill(0.0);
        for(fr=b;fr<=e;++fr) for(i=0; i<n; ++i) res.row(i) = system->traj[fr].coord[_index[i]];
    }
    res /= (e-b+1);
    return res;
}


////////////////////////////////////////////
// Transformations and inquery functions
////////////////////////////////////////////

bool Selection::is_large(){
    Vector3f _min,_max;
    minmax(_min,_max);
    return ( (box().lab_to_box(_max-_min) - 0.5*box().extents()).array() > 0).any();
}


void Selection::process_pbc_atom(int& a) const {
    if(a>=size()) throw PterosError("Wrong pbc atom {} for selection with {} atoms!",a,size());
    if(a<0) a = size()/2;
}


// Center of geometry
Vector3f Selection::center(bool mass_weighted, Array3i_const_ref pbc, int pbc_atom) const {    
    int n = size();

    if(n==0) throw PterosError("Can't get center of empty selection!");
    // in case of just one atom nothing to compute
    if(n==1) return xyz(0);

    Vector3f res;
    int i;

    process_pbc_atom(pbc_atom);

    res.fill(0.0);

    if( (pbc==0).all() ){
        // Non-periodic variant
        if(mass_weighted){
            float M = 0.0;
            #pragma omp parallel
            {                
                Vector3f r(Vector3f::Zero());
                #pragma omp for nowait reduction(+:M)
                for(i=0; i<n; ++i){
                    r += xyz(i)*mass(i);
                    M += mass(i);
                }
                #pragma omp critical
                {
                    res += r;
                }
            }
            if(M==0) throw PterosError("Selection has zero mass! Center of mass failed!");
            return res/M;
        } else {
            #pragma omp parallel
            {
                Vector3f r(Vector3f::Zero()); // local to omp thread
                #pragma omp for nowait
                for(i=0; i<n; ++i) r += xyz(i);
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
        // using leading point as a reference
        Vector3f ref_point = xyz(pbc_atom);
        PeriodicBox& b = system->box(frame);
        if(mass_weighted){
            float M = 0.0;
            #pragma omp parallel
            {                
                Vector3f r(Vector3f::Zero());
                #pragma omp for nowait reduction(+:M)
                for(i=0; i<n; ++i){
                    r += b.closest_image(xyz(i),ref_point,pbc) * mass(i);
                    M += mass(i);
                }
                #pragma omp critical
                {
                    res += r;                    
                }
            }
            if(M==0) throw PterosError("Selection has zero mass! Center of mass failed!");
            return res/M;
        } else {
            #pragma omp parallel
            {
                Vector3f r(Vector3f::Zero()); // local to omp thread
                #pragma omp for nowait
                for(i=0; i<n; ++i)
                    r += b.closest_image(xyz(i),ref_point,pbc);
                #pragma omp critical
                {
                    res += r;
                }
            }
            return res/n;
        }
    }
}

Vector3f Selection::center(const std::vector<float> &weights, Array3i_const_ref pbc, int pbc_atom) const
{
    int n = size();

    if(n==0) throw PterosError("Can't get center of empty selection!");
    if(weights.size()!=size()) throw PterosError("Weights array has {} elements, while selection has {}!",weights.size(),size());
    // in case of just one atom nothing to compute
    if(n==1) return xyz(0);

    Vector3f res;
    int i;

    process_pbc_atom(pbc_atom);

    res.fill(0.0);

    if( (pbc==0).all() ){
        // Non-periodic variant
        float W = 0.0;
        #pragma omp parallel
        {
            Vector3f r(Vector3f::Zero());
            #pragma omp for nowait reduction(+:W)
            for(i=0; i<n; ++i){
                r += xyz(i)*weights[i];
                W += weights[i];
            }
            #pragma omp critical
            {
                res += r;
            }
        }
        if(W==0) throw PterosError("Sum of weights is zero mass! Center failed!");
        return res/W;

    } else {
        // Periodic center
        // We will find closest periodic images of all points
        // using leading point as a reference
        Vector3f ref_point = xyz(pbc_atom);
        PeriodicBox& b = system->box(frame);

        float W = 0.0;
        #pragma omp parallel
        {
            Vector3f r(Vector3f::Zero());
            #pragma omp for nowait reduction(+:W)
            for(i=0; i<n; ++i){
                r += b.closest_image(xyz(i),ref_point,pbc) * weights[i];
                W += weights[i];
            }
            #pragma omp critical
            {
                res += r;
            }
        }
        if(W==0) throw PterosError("Sum of weights is zero mass! Center failed!");
        return res/W;
    }
}

// Plain translation
void Selection::translate(Vector3f_const_ref v){
    int i,n = _index.size();
    #pragma omp parallel for
    for(i=0; i<n; ++i) xyz(i) += v;
}

void Selection::translate_to(Vector3f_const_ref p, bool mass_weighted, Array3i_const_ref pbc, int pbc_atom){
    process_pbc_atom(pbc_atom);
    auto c = center(mass_weighted,pbc,pbc_atom);
    translate(p-c);
}

// Rotation around given vector relative to pivot
void Selection::rotate(Vector3f_const_ref pivot, Vector3f_const_ref axis, float angle){
    apply_transform( rotation_transform(pivot,axis,angle) );
}

////////////////////////////////////
// RMSD and fitting functions
////////////////////////////////////

// RMSD between two frames
float Selection::rmsd(int fr1, int fr2) const{
    int n = _index.size();
    float res = 0.0;

    if(fr1<0 || fr1>=system->num_frames() ||
       fr2<0 || fr2>=system->num_frames())
        throw PterosError("RMSD requested for frames {}:{} while the valid range is 0:{}", fr1,fr2,system->num_frames()-1);

    #pragma omp parallel for reduction(+:res)
    for(int i=0; i<n; ++i)
        res += (xyz(i,fr1)-xyz(i,fr2)).squaredNorm();

    return sqrt(res/n);
}

//RMSD between current and other frame
float Selection::rmsd(int fr) const {
    if(fr<0 || fr>=system->num_frames())
        throw PterosError("RMSD requested for frame {} while the valid range is 0:{}", fr,system->num_frames()-1);

    return rmsd(frame,fr);
}

// Apply transformation
void Selection::apply_transform(const Affine3f &t){
    int n = size();    
    #pragma omp parallel for
    for(int i=0; i<n; ++i){        
        xyz(i) = t * xyz(i);
    }
}

namespace pteros {

void copy_coord(const Selection &from, int from_fr, Selection &to, int to_fr)
{
    if(from.size()!=to.size())
        throw PterosError("Can't copy coordinates between selections of different size!");

    for(int i=0; i<from.size(); ++i){
        to.xyz(i,to_fr) = from.xyz(i,from_fr);
    }
}

void copy_coord(const Selection &from, Selection &to){
    copy_coord(from,from.get_frame(),to,to.get_frame());
}

Vector2f non_bond_energy(const Selection& sel1,
                         const Selection& sel2,
                         float cutoff,
                         int fr,
                         bool pbc)
{
    // Check if both selections are from the same system
    if(sel1.get_system()!=sel2.get_system())
        throw PterosError("Can't compute non-bond energy between selections from different systems!");

    if(fr<0) fr = sel1.get_frame();

    // Need to set frame fr for both selection to get correct grid search
    int fr1 = sel1.get_frame();
    int fr2 = sel2.get_frame();

    if(fr1!=fr) const_cast<Selection&>(sel1).set_frame(fr);
    if(fr2!=fr) const_cast<Selection&>(sel2).set_frame(fr);

    float d;
    if(cutoff==0){
        d = sel1.get_system()->get_force_field().get_cutoff();
    } else {
        d = cutoff;
    }

    // Perform grid search
    vector<Vector2i> pairs;
    vector<float> dist;
    Vector3i pbc_dims = pbc ? fullPBC : noPBC;
    search_contacts(d,sel1,sel2,pairs,dist,true, pbc_dims);

    // Restore frames if needed
    if(fr1!=fr) const_cast<Selection&>(sel1).set_frame(fr1);
    if(fr2!=fr) const_cast<Selection&>(sel2).set_frame(fr2);

    // Now get energy using pair list and distances
    return get_energy_for_list(pairs,dist,*sel1.get_system());
}


// RMSD between two selections (specified frames)
float rmsd(const Selection& sel1, int fr1, const Selection& sel2, int fr2){
    int n1 = sel1._index.size();
    int n2 = sel2._index.size();
    float res = 0.0;
    if(n1!=n2) throw PterosError("Incompatible selections for RMSD of sizes {} and {}", n1,n2);

    if(fr1<0 || fr1>=sel1.system->num_frames() || fr2<0 || fr2>=sel2.system->num_frames())
        throw PterosError("RMSD requested for frames {}:{} while the valid range is {}:{}",
                          fr1,fr2,sel1.system->num_frames()-1,sel2.system->num_frames()-1);


    #pragma omp parallel for reduction(+:res)
    for(int i=0; i<n1; ++i)
        res += (sel1.xyz(i,fr1)-sel2.xyz(i,fr2)).squaredNorm();

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

    if(n1!=n2) throw PterosError("Incompatible selections for fitting of sizes {} and {}", n1, n2);

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
            _u += sel1.xyz(i)*sel2.xyz(i).transpose()*sel1.mass(i);
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
    return Translation3f(cm2) * rot * Translation3f(-cm1);
}

// Fit two selection directly
void fit(Selection& sel1, const Selection& sel2){    
    Affine3f t = pteros::fit_transform(sel1,sel2);
    sel1.apply_transform(t);
}

} //namespace pteros


Vector2f Selection::non_bond_energy(float cutoff, bool pbc) const
{
    float d;
    if(cutoff==0){
        d = system->get_force_field().get_cutoff();
    } else {
        d = cutoff;
    }

    // Perform grid search
    vector<Vector2i> pairs;
    vector<float> dist;
    Vector3i pbc_dims = pbc ? fullPBC : noPBC;
    search_contacts(d,*this,pairs,dist,true, pbc_dims);

    // Now get energy using pair list and distances
    return get_energy_for_list(pairs,dist,*system);
}

// Fit all frames in trajectory
void Selection::fit_trajectory(int ref_frame, int b, int e){
    if(e==-1) e = system->num_frames()-1;
    // Sanity check
    if(b<0 || b>=system->num_frames() || b>e) throw PterosError("Invalid frame range!");
    if(ref_frame<0 || ref_frame>=system->num_frames())
        throw PterosError("Reference frame is out of range!");

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
    n = _index.size();

    float x_min, y_min, z_min, x_max, y_max, z_max;
    x_min = y_min = z_min = 1e10;
    x_max = y_max = z_max = -1e10;

    #pragma omp parallel for reduction(min:x_min,y_min,z_min) reduction(max:x_max,y_max,z_max)
    for(i=0; i<n; ++i){
        p = xyz_ptr(i);
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

void Selection::write(string fname, int b, int e) const {
    // -1 has special meaning
    if(b==-1) b=get_frame(); // current frame
    if(e==-1) e=system->num_frames()-1; // last frame

    int last = get_system()->num_frames();

    if(b<-1 || b>=last) throw PterosError("First frame {} for writing is out of range [{}:{}]",b,0,last);
    if(e<-1 || e>=last) throw PterosError("Last frame {} for writing is out of range [{}:{}]",e,0,last);
    if(e<b) throw PterosError("Last frame {} before first one {} for writing!",e,b);

    auto f = FileHandler::open(fname,'w');

    if(!(f->get_content_type().traj()) && e!=b){
        throw PterosError("Can't write the range of frames to structure file!");
    }    

    LOG()->debug("Writing {} frames to file '{}'...",e-b+1,fname);

    // Save current frame
    int cur_fr = get_frame();
    for(int fr=b;fr<=e;++fr){        
        const_cast<Selection*>(this)->set_frame(fr);
        f->write(*this,f->get_content_type());
    }
    const_cast<Selection*>(this)->set_frame(cur_fr);
}

void Selection::write(const std::unique_ptr<FileHandler> &handler, FileContent what, int b, int e) const
{
    // -1 has special meaning
    if(b==-1) b=get_frame(); // current frame
    if(e==-1) e=system->num_frames()-1; // last frame

    int last = get_system()->num_frames();

    if(b<-1 || b>=last) throw PterosError("First frame {} for writing is out of range [{}:{}]",b,0,last);
    if(e<-1 || e>=last) throw PterosError("Last frame {} for writing is out of range [{}:{}]",e,0,last);
    if(e<b) throw PterosError("Last frame {} before first one {} for writing!",e,b);

    if(!(handler->get_content_type().traj()) && e!=b && what.traj()){
        throw PterosError("Can't write the range of frames to this file!");
    }

    // First write all except trajectory (if any)
    if(what.atoms() || what.coord()){
        auto c = what;
        c.traj(false);
        handler->write(*this,c);
    }

    LOG()->debug("Writing {} frames to file handler...",e-b);

    // Now write trajectory if asked
    if(what.traj()){
        int cur_fr = get_frame();
        for(int fr=b;fr<=e;++fr){
            const_cast<Selection*>(this)->set_frame(fr);
            handler->write(*this, FileContent().traj(true));
        }
        const_cast<Selection*>(this)->set_frame(cur_fr);
    }
}

bool Selection::text_based() const {
    return sel_text!="";
}

bool Selection::coord_dependent() const {
    return (bool)parser;
}

void Selection::flatten()
{
    parser.reset();
    sel_text = "";
}

string Selection::to_gromacs_ndx(string name) const
{    
    stringstream s;
    s << "[ " << name << " ]" << endl;
    int n=0;
    for(int ind: _index){
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

int Selection::find_index(int global_index) const
{
    auto res = find(_index.begin(),_index.end(),global_index);
    if(res!=_index.end())
        return res-_index.begin();
    else
        return -1;
}


std::vector<std::vector<int>> Selection::get_internal_bonds(float d, bool periodic) const
{
    vector<vector<int>> con;

    if(d==0){
        // Use bonds from topology
        get_local_bonds_from_topology(con);        
    } else {
        // Find all connectivity pairs for given cut-off
        vector<Vector2i> pairs;
        vector<float> dist;
        Vector3i pbc_dims = periodic ? fullPBC : noPBC;
        search_contacts(d,*this,pairs,dist,false, pbc_dims); // local indexes

        // Form a connectivity structure in the form con[i]->1,2,5...
        con.resize(size());
        for(int i=0; i<pairs.size(); ++i){
            con[pairs[i](0)].push_back(pairs[i](1));
            con[pairs[i](1)].push_back(pairs[i](0));
        }
    }

    return con;
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
        ind = resindex(i);
        if(used.count(ind)) continue;
        // Starting global index
        b = e = index(i);
        // Go backward
        while( b-1>=0 && system->atoms[b-1].resindex == ind){ --b; };
        // Go forward
        while( e+1<system->atoms.size() && system->atoms[e+1].resindex == ind){ ++e; };
        sel.emplace_back(*system,b,e);
        // Mark as used
        used.insert(ind);
    }
}


void Selection::dssp(string fname) const {
    // Take all protein residues from current selection and select them whole
    auto sel = const_cast<Selection*>(this)->select("by residue protein");
    if(!sel.size()) PterosError("Need some protein residues for DSSP!");

    ofstream f(fname.c_str());
    dssp_wrapper(sel,f);
    f.close();
}

void Selection::dssp(ostream& os) const {
    // Take all protein residues from current selection and select them whole
    auto sel = const_cast<Selection*>(this)->select("by residue protein");
    if(!sel.size()) PterosError("Need some protein residues for DSSP!");

    dssp_wrapper(sel,os);
}


string Selection::dssp() const{
    // Take all protein residues from current selection and select them whole
    auto sel = const_cast<Selection*>(this)->select("by residue protein");
    if(!sel.size()) PterosError("Need some protein residues for DSSP!");

    return dssp_string(sel);
}


float Selection::vdw(int ind) const {
    return get_vdw_radius(atomic_number(ind),name(ind));
}

string Selection::element_name(int ind) const {
    int el = atomic_number(ind);
    return get_element_name(el);
}


MatrixXf Selection::atom_traj(int ind, int b, int e, bool make_row_major_matrix) const {
    if(e==-1) e = system->num_frames()-1;
    // Sanity check
    if(ind<0 || ind>=_index.size()) throw PterosError("Selection index is out of range!");
    if(b<0 || b>=system->num_frames() || b>e) throw PterosError("Invalid frame range!");

    int Nfr = e-b+1;
    MatrixXf ret;

    if(!make_row_major_matrix){
        ret.resize(3,Nfr);
        for(int fr=b;fr<=e;++fr) ret.col(fr) = system->traj[fr].coord[_index[ind]];
    } else {
        ret.resize(Nfr,3);
        for(int fr=b;fr<=e;++fr) ret.row(fr) = system->traj[fr].coord[_index[ind]];
    }

    return ret;
}

void Selection::split_by_connectivity(float d, std::vector<Selection> &res, bool periodic) {
    res.clear();

    vector<vector<int>> con = get_internal_bonds(d,periodic);

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
        tmp._index.push_back(_index[cur]);
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
                    res.back()._index.push_back(_index[i]);
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
                if(i>=size()) throw PterosError("UPS! Wrong connectivity?");
            }

            todo.insert(i);
            used[i] = 1;
            Nused++;
            // Start new selection
            Selection tmp(*system);
            tmp.set_frame(frame);
            tmp.sel_text = "";
            tmp._index.push_back(_index[i]);
            res.push_back(tmp);
        } else {
            // Done.
            // Indexes of selections are not sorted sort them
            for(auto& s: res) sort(s._index.begin(),s._index.end());
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
        m[resindex(i)].push_back(index(i));
    }
    // Create selections
    map<int,vector<int> >::iterator it;
    for(it=m.begin();it!=m.end();it++){
        res.emplace_back(*system, it->second.begin(), it->second.end());
    }
}

void Selection::split_by_molecule(std::vector<Selection> &res)
{
    if(!system->force_field.ready) throw PterosError("Can't split by molecule: no topology!");

    map<int,vector<int>> m;
    int bmol = 0;
    for(int i=0;i<size();++i){
        // Find molecule for this atom
        for(int j=bmol; j<system->force_field.molecules.size(); ++j){
            if(_index[i]>=system->force_field.molecules[j](0) && _index[i]<=system->force_field.molecules[j](1)){
                bmol=j;
                m[j].push_back(_index[i]);
            }
        }
    }

    // Create selections
    map<int,vector<int> >::iterator it;
    for(it=m.begin();it!=m.end();it++){
        res.emplace_back(*system, it->second.begin(), it->second.end());
    }
}

void Selection::split_by_chain(std::vector<Selection> &chains)
{
    // We split selection into several by chain
    chains.clear();
    // Map of resindexes to indexs in selections
    map<char,vector<int> > m;
    for(int i=0; i<size(); ++i){
        m[chain(i)].push_back(index(i));
    }
    // Create selections
    map<char,vector<int> >::iterator it;
    for(it=m.begin();it!=m.end();it++){
        chains.emplace_back(*system, it->second.begin(), it->second.end());
    }
}

void Selection::split_by_contiguous_index(std::vector<Selection> &parts)
{
    parts.clear();
    // Start first contiguous part
    int b = 0, i = 0;
    while(i<size()){
        while(i+1<size() && _index[i+1]==_index[i]+1) ++i;
        // Part finished
        parts.emplace_back(*system, _index[b],_index[i]);
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
        while(i+1<size() && (resindex(i+1)==resindex(i)+1 || resindex(i+1)==resindex(i)) ) ++i;
        // Part finished
        parts.emplace_back(*system,_index[b],_index[i]);
        b = i+1;
        i = b;
    }
}

void Selection::inertia(Vector3f_ref moments, Matrix3f_ref axes, Array3i_const_ref pbc, int pbc_atom) const{
    int n = size();
    int i;
    // Compute the central tensor of inertia. Place it into axes

    process_pbc_atom(pbc_atom);

    float axes00=0.0, axes11=0.0, axes22=0.0, axes01=0.0, axes02=0.0, axes12=0.0;

    Vector3f c = center(true,pbc,pbc_atom);

    if( (pbc!=0).any() ){
        Vector3f anchor = xyz(pbc_atom);
        PeriodicBox& b = system->box(frame);
        #pragma omp parallel
        {
            Vector3f p,d;
            float m;            
            #pragma omp for reduction(+:axes00,axes11,axes22,axes01,axes02,axes12)
            for(i=0;i<n;++i){
                // 0 point was used as an anchor in periodic center calculation,
                // so we have to use it as an anchor here as well!
                p = b.closest_image(xyz(i),anchor,pbc);
                d = p-c;
                m = mass(i);
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
                d = xyz(i)-c;
                m = mass(i);
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

float Selection::gyration(Array3i_const_ref pbc, int pbc_atom) const {
    process_pbc_atom(pbc_atom);

    int n = size();
    int i;
    float d, a = 0.0, b = 0.0;
    Vector3f c = center(true,pbc,pbc_atom);

    #pragma omp parallel for private(d) reduction(+:a,b)
    for(i=0;i<n;++i){
        if( (pbc!=0).any() ){
            d = system->box(frame).distance(xyz(i),c,pbc);
        } else {
            d = (xyz(i)-c).norm();
        }        
        a += mass(i)*d*d;
        b += mass(i);
    }
    return sqrt(a/b);
}

Vector3f Selection::dipole(bool as_charged, Array3i_const_ref pbc, int pbc_atom) const {
    process_pbc_atom(pbc_atom);

    Vector3f shift(0,0,0);
    if(as_charged) shift = center(true,pbc,pbc_atom);

    Vector3f res(0,0,0);
    for(int i=0;i<size();++i){
        if( (pbc!=0).any() ){
            res += (system->box(frame).closest_image(xyz(i),xyz(pbc_atom),pbc)-shift) * charge(i);
        } else {
            res += (xyz(i)-shift)*charge(i);
        }
    }
    return res * 0.02081943; // Convert to Debye
}

float Selection::distance(int i, int j, Array3i_const_ref pbc) const
{
    return system->distance(index(i),index(j),frame, pbc);
}

float Selection::angle(int i, int j, int k, Array3i_const_ref pbc) const
{
    return system->angle(index(i),index(j),index(k),frame,pbc);
}

float Selection::dihedral(int i, int j, int k, int l, Array3i_const_ref pbc) const
{
    return system->dihedral(index(i),index(j),index(k),index(l),frame,pbc);
}

void Selection::wrap(Array3i_const_ref pbc){
    for(int i=0;i<size();++i){
        box().wrap_point(xyz(i),pbc);
    }
}

void Selection::unwrap(Array3i_const_ref pbc, int pbc_atom){
    process_pbc_atom(pbc_atom);
    Vector3f c;
    if( (pbc==0).all() && pbc_atom>=0){
        c = xyz(pbc_atom);
    } else {
        c = center(true,pbc,pbc_atom);
    }
    for(int i=0;i<size();++i){
        xyz(i) = box().closest_image(xyz(i),c,pbc);
    }
}

int Selection::unwrap_bonds(float d, Array3i_const_ref pbc, int pbc_atom){
    process_pbc_atom(pbc_atom);
    int Nparts = 1;

    // a connectivity structure in the form con[i]->1,2,5...
    vector<vector<int>> con = get_internal_bonds(d,true); // periodic by definition

    // Mask of used atoms
    VectorXi used(size());
    used.fill(0);

    // Queue of centers to consider
    set<int> todo;

    // Leading atom coordinates
    Vector3f leading;

    todo.insert(pbc_atom);
    used[pbc_atom] = 1;
    int Nused = 1;
    PeriodicBox& b = system->box(frame);

    for(;;){
        while(!todo.empty()){
            // Get center from the queue
            int cur = *(todo.begin());
            todo.erase(todo.begin()); // And pop it from the queue
            leading = xyz(cur);

            // Find all atoms connected to this one
            for(int i=0; i<con[cur].size(); ++i){
                // We only add atoms, which were not yet used as centers
                if(used(con[cur][i])==0){
                    // Unwrap atom
                    xyz(con[cur][i]) = b.closest_image(xyz(con[cur][i]),leading,pbc);
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
                if(i>=size()) throw PterosError("Wrong connectivity?");
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

Eigen::Affine3f Selection::principal_transform(Array3i_const_ref pbc, int pbc_atom) const {
    process_pbc_atom(pbc_atom);
    Affine3f rot;
    Vector3f cm = center(true,pbc,pbc_atom);
    // We need to const-cast in order to call non-const translate() from const method
    // It's Ok because we'll translate back later
    const_cast<Selection*>(this)->translate(-cm);

    Matrix3f m;
    // Compute principal axes
    Vector3f moments;
    Matrix3f axes;
    inertia(moments,axes,pbc,pbc_atom);
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

void Selection::principal_orient(Array3i_const_ref pbc, int pbc_atom){
    Affine3f tr = principal_transform(pbc,pbc_atom);
    apply_transform(tr);
}

int Selection::num_residues() const
{
    if(size()==0) return 0;

    int count = 0;
    int cur = resindex(0);
    for(int i=1;i<size();++i){
        if(resindex(i)!=cur){
            cur = resindex(i);
            ++count;
        }
    }
    ++count; // Last one
    return count;
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
    reference operator*() const { return parent->xyz(pos); }
    pointer operator->() { return parent->xyz_ptr(pos); }
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
    for(int i=0; i<size(); ++i) radii[i] = vdw(i) + probe_r;

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
        //cout << i << " " << v << endl;
        if(do_a_per_atom) (*area_per_atom)[i] = v;
        surf += v;
    }

    return surf;
}

#else
float Selection::powersasa(float probe_r, vector<float> *area_per_atom,
                           float *total_volume, vector<float> *volume_per_atom) const
{
    throw PterosError("Pteros is compiled without powersasa support!");
}

#endif


// sasa implementation from MDTraj

float Selection::sasa(float probe_r, vector<float> *area_per_atom, int n_sphere_points) const
{
    // SASA code works for a range of frames. We use only current frame here

    // Make contigous array of coordinates
    vector<Vector3f> coord(size());
    for(int i=0;i<size();++i) coord[i] = xyz(i);

    vector<float> radii(size());
    for(int i=0; i<size(); ++i) radii[i] = vdw(i) + probe_r;

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


void Selection::get_local_bonds_from_topology(vector<vector<int>>& con) const {
    if(!system->force_field.ready) throw PterosError("No topology!");
    if(system->force_field.bonds.size()==0) throw PterosError("No bonds in topology!");

    con.clear();
    con.resize(size());

    int bind = index(0);
    int eind = index(size()-1);
    int a1,a2;
    auto bit = std::begin(_index);    
    auto eit = std::end(_index);

    for(int i=0;i<system->force_field.bonds.size();++i){
        a1 = system->force_field.bonds[i](0);
        a2 = system->force_field.bonds[i](1);
        if(a1>=bind && a1<=eind && a2>=bind && a2<=eind){
            auto it1 = std::find(bit,eit,a1);
            auto it2 = std::find(it1,eit,a2);
            con[it1-bit].push_back(it2-bit);
            con[it2-bit].push_back(it1-bit);
        }
    }
}



