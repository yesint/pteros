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



#include <fstream>
#include <unordered_map>
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "pteros/core/mol_file.h"
#include "selection_parser.h"
// DSSP
#include "pteros_dssp_wrapper.h"
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;
using namespace Eigen;
using namespace fmt;



// Base constructor of the system class
System::System() {

}

// Construnt system from file
System::System(string fname) {
    clear();
    load(fname);
}

System::System(const System& other){
    if(&other==this) return;
    clear();
    atoms = other.atoms;
    traj = other.traj;
    force_field = other.force_field;
}

System::System(const Selection &sel){
    if(sel.get_system()==this) throw Pteros_error("Can't construct system from selection of itself!");
    append(sel);
}


System& System::operator=(const System& other){
    if(&other==this) return *this;
    clear();
    atoms = other.atoms;
    traj = other.traj;
    force_field = other.force_field;
    return *this;
}

// Clear the system (i.e. before reading new system from file)
void System::clear(){
    atoms.clear();
    traj.clear();
    force_field.clear();
    filter.clear();
    filter_text = "";
}

void check_num_atoms_in_last_frame(const System& sys){
    if(sys.Frame_data(sys.num_frames()-1).coord.size()!=sys.num_atoms())
        throw Pteros_error("File contains {} atoms while the system has {}",
                       sys.Frame_data(sys.num_frames()-1).coord.size(), sys.num_atoms() );
}

void System::filter_atoms()
{
    if(filter.empty() && filter_text=="") return; // Do nothing if no filtering

    if(filter_text!="" && filter.empty()){
        // Parse text-based filter
        Selection_parser parser;
        parser.create_ast(filter_text);
        if(parser.has_coord) throw Pteros_error("Coordinate-dependent selections are not allowed in filters!");
        parser.apply(this, 0, filter);
    }

    vector<Atom> tmp = atoms;
    atoms.resize(filter.size());
    for(int i=0; i<filter.size(); ++i) atoms[i] = tmp[filter[i]];
}

void System::filter_coord(int fr)
{
    if(filter.empty()) return; // Do nothing if no filtering
    vector<Vector3f> tmp = traj[fr].coord;
    traj[fr].coord.resize(filter.size());
    for(int i=0; i<filter.size(); ++i) traj[fr].coord[i] = tmp[filter[i]];
}

// Load structure or trajectory
void System::load(string fname, int b, int e, int skip, std::function<bool(System*,int)> on_frame){
    // Create an IO handler
    auto f = Mol_file::open(fname,'r');

    int num_stored = 0;    

    if(num_atoms()==0){
        // We don't have atoms yet, so we will read everything possible except trajectory
        auto c = f->get_content_type();

        if(c.coord() && !c.traj()){
            // We have single frame. Read it directly here
            Frame fr;
            frame_append(fr);
            f->read(this, &Frame_data(num_frames()-1), c);

            filter_atoms();
            filter_coord(num_frames()-1);

            check_num_atoms_in_last_frame(*this);
            ++num_stored;

            assign_resindex();

            // Call a callback if asked
            if(on_frame) on_frame(this,num_frames()-1);            
            // And now we should just exit
            return;
        } else {
            // We have not a single frame, so read only atoms here
            // Clear flags for trajectory and coordinates
            c.coord(false);
            c.traj(false);

            f->read(this, nullptr, c);
            filter_atoms();
            assign_resindex();
        }        
    }

    if(num_atoms()>0){
        // We have atoms already, so read only coordinates
        auto c = f->get_content_type();
        if(!c.coord() && !c.traj())
            throw Pteros_error("File reader for file '{}' is not capable of appending frames to the system!",fname);

        // Check if we can read trajectory
        if(c.traj()){
            // Sanity check for frame range
            if((e<b && e!=-1)|| b<0)
                throw Pteros_error("Invalid frame range for reading!");

            int cur = 0; // This holds real frame index in trajectory

            // Skip frames if needed
            if(b>0){
                LOG()->info("Skipping {} frames...\n", b);
                Frame skip_fr;
                for(int i=0;i<b;++i){
                    f->read(nullptr, &skip_fr, Mol_file_content().traj(true));
                    cur++;
                }
            }            

            LOG()->info("Reading...\n");

            int actually_read = 0;

            bool callback_ok = true;

            while(true){
                // End frame reached?
                if(cur==e && e!=-1) break;

                // Append new frame where the data will go
                Frame fr;
                frame_append(fr);
                // Try to read into it
                bool ok = f->read(nullptr, &Frame_data(num_frames()-1), Mol_file_content().traj(true));
                if(!ok){
                    frame_delete(num_frames()-1); // Remove last frame - it's invalid
                    break; // end of trajectory
                }

                filter_coord(num_frames()-1);
                check_num_atoms_in_last_frame(*this);

                ++cur;
                ++actually_read;

                if(skip>0 && actually_read%skip!=0){
                    frame_delete(num_frames()-1); // Remove last frame - it's invalid
                    continue;
                } else {
                    actually_read = 0;
                }

                // If we are here new frame is accepted
                ++num_stored;

                // Call a callback if asked
                if(on_frame) callback_ok = on_frame(this,num_frames()-1);
                // If callback returns false stop reading
                if(!callback_ok) break;
            }
        } else if(c.coord() && !c.top()) {
            // File contains single frame and no topology
            // Append new frame where the data will go
            Frame fr;
            frame_append(fr);
            // Read it            
            f->read(nullptr, &Frame_data(num_frames()-1), Mol_file_content().coord(true));
            filter_coord(num_frames()-1);
            check_num_atoms_in_last_frame(*this);
            ++num_stored;
            // Call a callback if asked
            if(on_frame) on_frame(this,num_frames()-1);

        } else if(f->get_content_type().top()) {
            // This is topology file, read only topology
            f->read(this, nullptr, Mol_file_content().top(true) );
            // For topology filtering is not possible
            if(!filter.empty() || filter_text!="") throw Pteros_error("Filtering is not possible when reading topology!");
        }
    }

    LOG()->info("Accepted {} frames. Now {} frames in the System.\n", num_stored, num_frames());
}

bool System::load(const std::unique_ptr<Mol_file>& handler, Mol_file_content what, std::function<bool (System *, int)> on_frame)
{        
    // Asked for structure or topology
    if(what.atoms() || what.top() || what.coord()){
        if(what.top() && (!filter.empty() || filter_text!=""))
            throw Pteros_error("Filtering is not possible when reading topology!");

        // We don't want to read traj here, so disable it even if asked to read it
        auto c = what;
        c.traj(false);

        // If we are reading new atoms we have to clear the system first
        if(what.atoms()) clear();

        // See if we asked for coordinates
        if(what.coord()){
            Frame fr;
            frame_append(fr);
            handler->read(this, &Frame_data(num_frames()-1), c);

            filter_atoms();
            filter_coord(num_frames()-1);

            check_num_atoms_in_last_frame(*this);
            if(what.atoms()) assign_resindex();
            // Call a callback if asked
            if(on_frame) on_frame(this,num_frames()-1);
        } else {
            // Not asked for coordinates
            handler->read(this, nullptr, c);
            filter_coord(num_frames()-1);
            if(what.atoms()) assign_resindex();
        }
    }

    // Asked for trajectory frame
    if(what.traj()){
        // Append new frame where the data will go
        Frame fr;
        frame_append(fr);
        // Try to read into it
        bool ok = handler->read(nullptr, &Frame_data(num_frames()-1), Mol_file_content().traj(true));
        if(!ok){
            frame_delete(num_frames()-1); // Remove last frame - it's invalid
            return false;
        }

        filter_coord(num_frames()-1);
        check_num_atoms_in_last_frame(*this);

        // Call a callback if asked
        if(on_frame) on_frame(this,num_frames()-1);
        // Since we are reading one frame return doesn't matter
    }

    return true;
}

std::vector<std::pair<string,Selection>> System::load_gromacs_ndx(string fname)
{
    stringstream ss;
    string dum, grname="", line;
    int ind;
    std::vector<std::pair<string,Selection>> res;
    vector<int> vec;

    ifstream f(fname);
    if(!f) throw Pteros_error("Can't open ndx file '{}'!",fname);
    while(getline(f,line)){
        ss.clear(); ss.str(line);
        // Check if this is group name
        if(line.find_first_of('[')!=std::string::npos){
            if(grname!="") res.back().second.modify(vec);
            ss >> dum >> grname;
            res.push_back({grname,Selection(*this)});
            vec.clear();
        } else {
            // not a name. Try getting numbers
            for(;;) {
                 ss >> ind;
                 if(ss.good())
                    vec.push_back(ind-1);
                 else
                     break;
            }
        }
    }
    if(grname!="") res.back().second.modify(vec);
    if(res.empty()) throw Pteros_error("File '{}' contains no index groups!",fname);
    return res;
}

void System::set_filter(string str)
{
    if(atoms.size()) throw Pteros_error("Filter could be set to empty system only!");
    filter_text = str;
    filter.clear();
}

void sort_and_remove_duplicates(std::vector<int>& index)
{
    if(index.size()){
        sort(index.begin(),index.end());
        vector<int>::iterator it = unique(index.begin(), index.end());
        index.resize( it - index.begin() );
        if(index[0]<0) throw Pteros_error("Negative index {} present in filter!",index[0]);
    }
}


void System::set_filter(int ind1, int ind2)
{
    if(atoms.size()) throw Pteros_error("Filter could be set to empty system only!");
    filter.reserve(ind2-ind1+1);
    for(int i=ind1; i<=ind2; ++i) filter.push_back(i);
    sort_and_remove_duplicates(filter);
    filter_text="";
}

void System::set_filter(const std::vector<int> &ind)
{
    if(atoms.size()) throw Pteros_error("Filter could be set to empty system only!");
    filter = ind;
    sort_and_remove_duplicates(filter);
    filter_text="";
}

void System::set_filter(vector<int>::iterator it1, vector<int>::iterator it2)
{
    if(atoms.size()) throw Pteros_error("Filter could be set to empty system only!");
    copy(it1,it2,back_inserter(filter));
    sort_and_remove_duplicates(filter);
    filter_text="";
}

// Destructor of the system class
System::~System() {}

int System::frame_dup(int fr){
    if(fr<0 || fr>=traj.size())
        throw Pteros_error("Invalid frame {} for duplication!",fr);
    traj.push_back(traj[fr]);
    return traj.size()-1;
}

void System::frame_copy(int fr1, int fr2){
    if(fr1<0 || fr1>=traj.size() || fr2<0 || fr2>=traj.size())
        throw Pteros_error("Invalid frames {} and {} for copying!",fr1,fr2);
    traj[fr2] = traj[fr1];    
}

// Delete the range of frames. e = -1 is default
void System::frame_delete(int b, int e){
    int i;    

    if(e==-1) e = num_frames()-1;
    if(e<b || b<0 || e>num_frames()-1) throw Pteros_error("Invalid frame range {}:{} for deletion",b,e);

    // Get iterators for deletion
    vector<Frame>::iterator b_it, e_it;
    b_it = traj.begin();
    for(i=0;i<b;++i) b_it++;
    e_it = b_it;
    for(;i<e;++i) e_it++;
    e_it++; //Go one past the end

    traj.erase(b_it,e_it);

    // Check if there are some frames left. If not print the warning
    // that all selections are invalid!
    if(traj.size()==0) LOG()->warn("All frames are deleted. All selections are now INVALID!\n");
}

void System::frame_swap(int fr1, int fr2)
{
    if(fr1<0 || fr1>=traj.size() || fr2<0 || fr2>=traj.size())
        throw Pteros_error("Invalid frames {} and {} for swapping!",fr1,fr2);
    std::swap(traj[fr1],traj[fr2]);
}

void System::frame_append(const Frame& fr){
    traj.push_back(fr);
}

void System::assign_resindex(int start){
    if(start<0) start=0;

    int curres = atoms[start].resid;
    int curchain = atoms[start].chain;
    int cur = 0;
    if(start>0) cur = atoms[start].resindex;
    for(int i=start; i<atoms.size(); ++i){
        if( atoms[i].resid!=curres || atoms[i].chain!=curchain ){
            ++cur;
            curres = atoms[i].resid;
            curchain = atoms[i].chain;
        }
        atoms[i].resindex = cur;
    }
}

void System::sort_by_resindex()
{
    // Make and array of indexes to shuffle
    vector<int> ind(atoms.size());
    for(int i=0;i<ind.size();++i) ind[i] = i;
    // Sort indexes
    sort(ind.begin(),ind.end(),
       [this](int i, int j){
        if(Atom_data(i).resindex == Atom_data(j).resindex){
            return (i<j);
        } else {
            return (Atom_data(i).resindex < Atom_data(j).resindex);
        }
       }
    );
    // Now shuffle atoms and coordinates according to indexes
    vector<Atom> tmp(atoms); //temporary
    for(int i=0;i<ind.size();++i) atoms[i] = tmp[ind[i]];

    std::vector<Eigen::Vector3f> tmp_coord;
    for(int j=0; j<traj.size(); ++j){ // Over all frames
        tmp_coord = traj[j].coord; //temporary
        for(int i=0;i<ind.size();++i) traj[j].coord[i] = tmp_coord[ind[i]];
    }
}

Selection System::atoms_dup(const vector<int>& ind){
    // Sanity check
    if(!ind.size()) throw Pteros_error("No atoms to duplicate!");
    for(int i=0; i<ind.size(); ++i){
        if(ind[i]<0 || ind[i]>atoms.size()-1)
            throw Pteros_error("Invalid index {} for atom duplication!", ind[i]);
    }

    // Duplicate atoms
    int first_added = atoms.size();
    int last_added = atoms.size()+ind.size()-1;
    // Prepare by increasing capacity of vectors
    atoms.reserve(atoms.size()+ind.size());
    for(int j=0; j<traj.size(); ++j){
        traj[j].coord.reserve(atoms.size()+ind.size());
    }

    // Now add atoms
    for(int i=0; i<ind.size(); ++i){
        // Add new atom
        atoms.push_back(atoms[ind[i]]);
        // Add new coordinate slot
        for(int j=0; j<traj.size(); ++j){
            traj[j].coord.push_back(traj[j].coord[ind[i]]);
        }
    }


    return Selection(*this,first_added,last_added);
}

Selection System::atoms_add(const vector<Atom>& atm, const vector<Vector3f>& crd){
    // Sanity check
    if(!atm.size()) throw Pteros_error("No atoms to add!");
    if(atm.size()!=crd.size())
        throw Pteros_error("Number of coordinates {} doesn't mutch number of atoms {}",crd.size(),atm.size());

    int first_added = atoms.size();
    int last_added = atoms.size()+atm.size()-1;
    // Prepare by increasing capacity of vectors
    atoms.reserve(atoms.size()+atm.size());

    // If no frames add one
    if(traj.size()==0) traj.push_back(Frame());

    for(int j=0; j<traj.size(); ++j){
        traj[j].coord.reserve(atoms.size()+atm.size());
    }
    // Now add atoms
    for(int i=0; i<atm.size(); ++i){
        // Add new atom
        atoms.push_back(atm[i]);
        // Add new coordinate slot
        for(int j=0; j<traj.size(); ++j){
            traj[j].coord.push_back(crd[i]);
        }
    }

    return Selection(*this,first_added,last_added);
}

void System::atoms_delete(const std::vector<int> &ind){
    int i,fr;

    // Sanity check
    if(!ind.size()) throw Pteros_error("No atoms to delete!");
    for(int i=0; i<ind.size(); ++i){
        if(ind[i]<0 || ind[i]>atoms.size()-1)
            throw Pteros_error("Index {} for atom is out of range (0:{})!",ind[i],atoms.size()-1);
    }

    // Mark atoms for deletion by assigning negative mass
    for(i=0;i<ind.size();++i)
        atoms[ind[i]].mass = -1.0;

    // Cycle over all atoms and keep only those with positive mass
    vector<pteros::Atom> tmp = atoms;
    atoms.clear();
    for(i=0;i<tmp.size();++i){
        if(tmp[i].mass>=0) atoms.push_back(tmp[i]);
    }

    // Now cycle over trajectory and keep only corresponding coordinates
    vector<Vector3f> tmp_coord;
    for(fr=0; fr<num_frames(); ++fr){
        // Make a copy of traj coords
        tmp_coord = traj[fr].coord;
        traj[fr].coord.clear();
        for(i=0;i<tmp.size();++i){
            if(tmp[i].mass>=0) traj[fr].coord.push_back(tmp_coord[i]);
        }
    }
}

void System::atom_move(int i, int j)
{
    // Sanity check
    if(i<0 || i>=num_atoms()) throw Pteros_error(format("Index of atom to move ({}}) is out of range ({}:{})!", i,0,num_atoms()));
    if(j<0 || i>=num_atoms()) throw Pteros_error(format("Target index to move ({}}) is out of range ({}:{}})!", i,0,num_atoms()));
    if(i==j) return; // Nothing to do

    // Move atom
    auto at = Atom_data(i);

    if(i<j){
        for(int a=i+1; a<=j; ++a) Atom_data(a-1) = Atom_data(a);
        Atom_data(j) = at;

        for(int fr=0; fr<num_frames(); ++fr){
            auto coord = XYZ(i,fr);
            for(int a=i+1; a<=j; ++a) XYZ(a-1,fr) = XYZ(a,fr);
            XYZ(j,fr) = coord;
        }
    } else {
        for(int a=i-1; a>=j; --a) Atom_data(a+1) = Atom_data(a);
        Atom_data(j) = at;

        for(int fr=0; fr<num_frames(); ++fr){
            auto coord = XYZ(i,fr);
            for(int a=i-1; a>=j; --a) XYZ(a+1,fr) = XYZ(a,fr);
            XYZ(j,fr) = coord;
        }
    }
}

Selection System::atom_clone(int source)
{
    atoms_dup({source});
    atom_move(num_atoms()-1,source+1);
    return Selection(*this,source+1,source+1);
}

Selection System::append(const System &sys){
    //Sanity check
    if(num_frames()>0 && num_frames()!=sys.num_frames())
        throw Pteros_error("Can't merge systems with different number of frames ({} and {}})!",num_frames(),sys.num_frames());
    // If no frames create needed ammount
    bool transfer_time_box = false;
    if(num_frames()==0){
        transfer_time_box = true;
        traj.resize(sys.num_frames());
    }

    int first_added = num_atoms();

    // Merge atoms
    copy(sys.atoms.begin(),sys.atoms.end(),back_inserter(atoms));
    // Merge coordinates
    for(int fr=0; fr<num_frames(); ++fr){
        if(transfer_time_box){
            traj[fr].time = sys.Time(fr);
            traj[fr].box = sys.Box(fr);
        }
        copy(sys.traj[fr].coord.begin(),sys.traj[fr].coord.end(),back_inserter(traj[fr].coord));
    }

    // Reassign resindex for added atoms
    assign_resindex(first_added-1);

    return Selection(*this,first_added,num_atoms()-1);
}

Selection System::append(const Selection &sel){    
    //Sanity check    
    if(num_frames()>0 && num_frames()!=sel.get_system()->num_frames())
        throw Pteros_error("Can't merge systems with different number of frames! ({} and {}})",
                           num_frames(),sel.get_system()->num_frames());

    // If no frames create needed ammount
    bool transfer_time_box = false;
    if(num_frames()==0){
        transfer_time_box = true;
        traj.resize(sel.get_system()->num_frames());        
    }

    int first_added = num_atoms();    

    // Merge atoms
    atoms.reserve(atoms.size()+sel.size());
    for(int i=0;i<sel.size();++i) atoms.push_back(sel.Atom_data(i));
    // Merge coordinates    
    for(int fr=0; fr<num_frames(); ++fr){
        traj[fr].coord.reserve(atoms.size());
        if(transfer_time_box){
            traj[fr].time = sel.get_system()->Time(fr);
            traj[fr].box = sel.get_system()->Box(fr);
        }
        for(int i=0;i<sel.size();++i){
            traj[fr].coord.push_back(sel.XYZ(i,fr));
        }
    }
    // Reassign resindex
    assign_resindex(first_added-1);

    return Selection(*this,first_added,num_atoms()-1);
}

Selection System::append(const Atom &at, Vector3f_const_ref coord)
{
    // If no frames create one
    if(num_frames()==0){
        traj.resize(1);
    }

    int first_added = num_atoms();

    // Add atom
    atoms.push_back(at);

    // Merge coordinates
    for(int fr=0; fr<num_frames(); ++fr){
        traj[fr].coord.push_back(coord);
    }
    // Reassign resindex
    assign_resindex(first_added-1);

    return Selection(*this,first_added,num_atoms()-1);
}

Selection System::append(const Atom_proxy &at)
{
    return append(at.atom(),at.xyz());
}

void System::rearrange(const std::vector<string> &sel_strings){
    vector<Selection> sel_vec(sel_strings.size());
    for (int i=0; i<sel_strings.size(); ++i){
        sel_vec[i].modify(*this, sel_strings[i]);        
    }

    rearrange(sel_vec);
}

void System::rearrange(const std::vector<Selection> &sel_vec){
    // Sanity check
    for(auto &s: sel_vec){
        if(s.size()==0) LOG()->warn("Empty selection in rearrange will be ignored!");
        if(s.get_system()!=this) throw Pteros_error("Rearrange needs selections from the same system!");
    }
    // Overlap check
    vector<int> inters;
    for (int i=0; i<sel_vec.size()-1; ++i){
        for (int j=i+1; j<sel_vec.size(); ++j){
            if(sel_vec[i].size()>0 && sel_vec[j].size()>0){
                set_intersection(sel_vec[i].index_begin(), sel_vec[i].index_end(),
                                 sel_vec[j].index_begin(), sel_vec[j].index_end(),
                                 back_inserter(inters));
            }
            if (!inters.empty()) throw Pteros_error("Selections for rearrange should not overlap!");
        }
    }

    Selection rest(*this);
    System result;

    // Append all explicitly given selections to result
    for (auto &s: sel_vec){
        if(s.size()>0){
            result.append(s);
            rest.append(s);
        }
    }
    // Invert to get rest
    rest = ~rest;

    // Only append rest to result if it is not empty
    if(rest.size()) result.append(rest);

    *this = result;
}

void System::keep(const string &sel_str)
{
    Selection sel(*this,sel_str);
    keep(sel);
}

void System::keep(const Selection &sel)
{
    if(sel.get_system()!=this) throw Pteros_error("keep needs selection from the same system!");
    System tmp;
    tmp.append(sel);
    *this = tmp;
}

void System::remove(const string &sel_str)
{    
    Selection sel(*this,"not ("+sel_str+")");
    keep(sel);
}

void System::remove(Selection &sel)
{
    if(sel.get_system()!=this) throw Pteros_error("remove needs selection from the same system!");
    System tmp;
    tmp.append(~sel);
    *this = tmp;
    sel.clear();
}


void System::distribute(const Selection sel, Vector3i_const_ref ncopies, Matrix3f_const_ref shift)
{
    if(sel.get_system()!=this) throw Pteros_error("distribute needs selection from the same system!");
    Selection tmp;
    Vector3f v;
    int n = sel.Resindex(sel.size()-1)+1;
    for(int x=0; x<ncopies(0); ++x)
        for(int y=0; y<ncopies(1); ++y)
            for(int z=0; z<ncopies(2); ++z){
                if(x>0 || y>0 || z>0){
                    tmp = append(sel);
                    v = shift.col(0)*x + shift.col(1)*y + shift.col(2)*z;
                    tmp.translate(v);
                    // Increment all resindexes of added selection
                    for(int i=0;i<tmp.size();++i) tmp.Resindex(i) += n;                    
                }
            }
}


Selection System::select(string str, int fr){
    return Selection(*this,str,fr);
}

Selection System::operator()(string str, int fr)
{
    return Selection(*this,str,fr);
}

Selection System::select(int ind1, int ind2){
    return Selection(*this,ind1,ind2);
}

Selection System::operator()(int ind1, int ind2)
{
    return Selection(*this,ind1,ind2);
}

Selection System::select(const std::vector<int> &ind){
    return Selection(*this,ind);
}

Selection System::operator()(const std::vector<int> &ind)
{
    return Selection(*this,ind);
}

Selection System::select(std::vector<int>::iterator it1, std::vector<int>::iterator it2){
    return Selection(*this,it1,it2);
}

Selection System::operator()(vector<int>::iterator it1, vector<int>::iterator it2)
{
    return Selection(*this,it1,it2);
}

Selection System::select(const std::function<void (const System &, int, std::vector<int> &)> &callback, int fr)
{
    return Selection(*this,callback,fr);
}

Selection System::operator()(const std::function<void (const System &, int, std::vector<int> &)> &callback, int fr)
{
    return Selection(*this,callback,fr);
}

Selection System::select_all(){
    return Selection(*this,0,num_atoms()-1);
}

Selection System::operator()()
{
    return select_all();
}

inline void wrap_coord(Vector3f& point, const Matrix3f& box,
                       const Vector3i dims_to_wrap = Vector3i::Ones()){
    Matrix3f b;
    b.col(0) = box.col(0).normalized();
    b.col(1) = box.col(1).normalized();
    b.col(2) = box.col(2).normalized();

    int i;
    float intp,fracp;
    // Get scalar projections onto box basis vectors
    Vector3f prj;    
    Vector3f box_dim = box.colwise().norm();

    prj = b.inverse()*point;

    for(i=0;i<3;++i){
        if(dims_to_wrap(i)!=0){
            fracp = std::modf(prj(i)/box_dim(i),&intp);
            if(fracp<0) fracp = fracp+1;
            prj(i) = box_dim(i)*fracp;
        }
    }   

    // Get back to lab coordinates
    point = b*prj;
}

float System::distance(int i, int j, int fr, Array3i_const_ref pbc) const {
    if( (pbc!=0).any() ){
        return traj[fr].box.distance(traj[fr].coord[i], traj[fr].coord[j], pbc);
    } else {
        return (traj[fr].coord[i] - traj[fr].coord[j]).norm();
    }
}


float System::angle(int i, int j, int k, int fr, Array3i_const_ref pbc) const
{
    Vector3f v1,v2;
    if( (pbc!=0).any() ){
        v1 = Box(fr).shortest_vector(XYZ(i,fr),XYZ(j,fr),pbc);
        v2 = Box(fr).shortest_vector(XYZ(k,fr),XYZ(j,fr),pbc);
    } else {
        v1 = XYZ(i,fr)-XYZ(j,fr);
        v2 = XYZ(k,fr)-XYZ(j,fr);
    }
    return acos(v1.dot(v2)/(v1.norm()*v2.norm()));
}

float System::dihedral(int i, int j, int k, int l, int fr, Array3i_const_ref pbc) const
{
    Vector3f b1,b2,b3;
    if( (pbc!=0).any() ){
        Vector3f _i = XYZ(i,fr);
        Vector3f _j = Box(fr).closest_image(XYZ(j,fr),_i,pbc);
        Vector3f _k = Box(fr).closest_image(XYZ(k,fr),_i,pbc);
        Vector3f _l = Box(fr).closest_image(XYZ(l,fr),_i,pbc);
        b1 = _j - _i;
        b2 = _k - _j;
        b3 = _l - _k;
    } else {
        b1 = XYZ(j,fr)-XYZ(i,fr);
        b2 = XYZ(k,fr)-XYZ(j,fr);
        b3 = XYZ(l,fr)-XYZ(k,fr);
    }    

    // Dihedral
    return atan2( ((b1.cross(b2)).cross(b2.cross(b3))).dot(b2/b2.norm()) ,
                  (b1.cross(b2)).dot(b2.cross(b3)) );
}

void System::wrap(int fr, Array3i_const_ref pbc){
    for(int i=0;i<num_atoms();++i){
        traj[fr].box.wrap_point(XYZ(i,fr),pbc);
    }
}

inline float LJ_en_kernel(float C6, float C12, float r_inv){
    //float tmp = 1/r;
    float tmp = r_inv*r_inv; // (1/r)^2
    tmp = tmp*tmp*tmp; // (1/r)^6
    return C12*tmp*tmp-C6*tmp;
}

#define ONE_4PI_EPS0      138.935456

inline float Coulomb_en_kernel(float q1, float q2, float r_inv){
    return ONE_4PI_EPS0*q1*q2*r_inv;
}

void System::dssp(string fname, int fr) const {
    ofstream f(fname.c_str());
    Selection sel(const_cast<System&>(*this),std::string("all"),fr);
    dssp_wrapper(sel,f);
    f.close();
}

void System::dssp(ostream& os, int fr) const {
    Selection sel(const_cast<System&>(*this),std::string("all"),fr);
    dssp_wrapper(sel,os);
}


string System::dssp(int fr) const{
    Selection sel(const_cast<System&>(*this),std::string("all"),fr);
    return dssp_string(sel);
}

