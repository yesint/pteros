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
#include <unordered_map>
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "pteros/core/file_handler.h"
#include "pteros/core/utilities.h"
#include "selection_parser.h"
#include "pteros/core/logging.h"
#include <utility>

using namespace std;
using namespace pteros;
using namespace Eigen;
using namespace fmt;

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

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
    if(sel.get_system()==this) throw PterosError("Can't construct system from selection of itself!");
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
    if(sys.frame(sys.num_frames()-1).coord.size()!=sys.num_atoms())
        throw PterosError("File contains {} atoms while the system has {}",
                       sys.frame(sys.num_frames()-1).coord.size(), sys.num_atoms() );
}

void System::filter_atoms()
{
    if(filter.empty() && filter_text=="") return; // Do nothing if no filtering

    if(filter_text!="" && filter.empty()){
        // Parse text-based filter
        SelectionParser parser;
        parser.create_ast(filter_text,this);
        if(parser.has_coord) throw PterosError("Coordinate-dependent selections are not allowed in filters!");
        parser.apply_ast(0, filter);
    }

    vector<Atom> tmp = atoms;
    atoms.resize(filter.size());
    for(int i=0; i<filter.size(); ++i) atoms[i] = tmp[filter[i]];
}

void System::filter_coord(int fr)
{
    if(filter.empty()) return; // Do nothing if no filtering
    Frame tmp = traj[fr];

    traj[fr].coord.resize(filter.size());    
    for(int i=0; i<filter.size(); ++i) traj[fr].coord[i] = tmp.coord[filter[i]];

    if(traj[fr].has_vel()){
        traj[fr].vel.resize(filter.size());
        for(int i=0; i<filter.size(); ++i) traj[fr].vel[i] = tmp.vel[filter[i]];
    }
    if(traj[fr].has_force()){
        traj[fr].force.resize(filter.size());
        for(int i=0; i<filter.size(); ++i) traj[fr].force[i] = tmp.force[filter[i]];
    }
}

// Load structure or trajectory
void System::load(string fname, int b, int e, int skip, std::function<bool(System*,int)> on_frame){
    // Create an IO handler
    auto f = FileHandler::open(fname,'r');

    int num_stored = 0;    

    if(num_atoms()==0){
        // We don't have atoms yet, so we will read everything possible except trajectory
        auto c = f->get_content_type();

        if(c.coord() && !c.traj()){
            // We have single frame. Read it directly here
            Frame fr;
            frame_append(fr);
            f->read(this, &frame(num_frames()-1), c);

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
            throw PterosError("File reader for file '{}' is not capable of appending frames to the system!",fname);

        // Check if we can read trajectory
        if(c.traj()){
            // Sanity check for frame range
            if((e<b && e!=-1)|| b<0)
                throw PterosError("Invalid frame range for reading!");

            int cur = 0; // This holds real frame index in trajectory

            // Skip frames if needed
            if(b>0){
                LOG()->info("Skipping {} frames...", b);
                Frame skip_fr;
                for(int i=0;i<b;++i){
                    f->read(nullptr, &skip_fr, FileContent().traj(true));
                    cur++;
                }
            }            

            LOG()->debug("Reading trajectory '{}'...",fname);

            int actually_read = 0;

            bool callback_ok = true;

            while(true){
                // End frame reached?
                if(cur==e && e!=-1) break;

                // Append new frame where the data will go
                Frame fr;
                frame_append(fr);
                // Try to read into it
                bool ok = f->read(nullptr, &frame(num_frames()-1), FileContent().traj(true));                

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
            f->read(nullptr, &frame(num_frames()-1), FileContent().coord(true));
            filter_coord(num_frames()-1);
            check_num_atoms_in_last_frame(*this);
            ++num_stored;
            // Call a callback if asked
            if(on_frame) on_frame(this,num_frames()-1);

        } else if(f->get_content_type().top()) {
            // This is topology file, read only topology
            f->read(this, nullptr, FileContent().top(true) );
            // For topology filtering is not possible
            if(!filter.empty() || filter_text!="") throw PterosError("Filtering is not possible when reading topology!");
        }
    }

    LOG()->info("Accepted {} frames. Now {} frames in the System.", num_stored, num_frames());
}

bool System::load(const std::unique_ptr<FileHandler>& handler, FileContent what, std::function<bool (System *, int)> on_frame)
{        
    // Asked for structure or topology
    if(what.atoms() || what.top() || what.coord()){
        if(what.top() && (!filter.empty() || filter_text!=""))
            throw PterosError("Filtering is not possible when reading topology!");

        // We don't want to read traj here, so disable it even if asked to read it
        auto c = what;
        c.traj(false);

        // If we are reading new atoms we have to clear the system first
        if(what.atoms()) clear();

        // See if we asked for coordinates
        if(what.coord()){
            Frame fr;
            frame_append(fr);
            handler->read(this, &frame(num_frames()-1), c);

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
        bool ok = handler->read(nullptr, &frame(num_frames()-1), FileContent().traj(true));
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

void System::write(string fname, int b, int e) const
{
    select_all().write(fname,b,e);
}

void System::write(const std::unique_ptr<FileHandler> &handler, FileContent what, int b, int e) const
{
    select_all().write(handler,what,b,e);
}

std::vector<std::pair<string,Selection>> System::load_gromacs_ndx(string fname)
{
    stringstream ss;
    string dum, grname="", line;
    int ind;
    std::vector<std::pair<string,Selection>> res;
    vector<int> vec;

    ifstream f(fname);
    if(!f) throw PterosError("Can't open ndx file '{}'!",fname);
    while(getline(f,line)){
        ss.clear(); ss.str(line);
        // Check if this is group name
        if(line.find_first_of('[')!=std::string::npos){
            if(grname!="") res.back().second.modify(vec);
            ss >> dum >> grname;
            res.emplace_back(grname,Selection(*this));
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
    if(res.empty()) throw PterosError("File '{}' contains no index groups!",fname);
    return res;
}

void System::set_filter(string str)
{
    if(atoms.size()) throw PterosError("Filter could be set to empty system only!");
    filter_text = str;
    filter.clear();
}

void sort_and_remove_duplicates(std::vector<int>& index)
{
    if(index.size()){
        sort(index.begin(),index.end());
        vector<int>::iterator it = unique(index.begin(), index.end());
        index.resize( it - index.begin() );
        if(index[0]<0) throw PterosError("Negative index {} present in filter!",index[0]);
    }
}


void System::set_filter(int ind1, int ind2)
{
    if(atoms.size()) throw PterosError("Filter could be set to empty system only!");
    filter.reserve(ind2-ind1+1);
    for(int i=ind1; i<=ind2; ++i) filter.push_back(i);
    sort_and_remove_duplicates(filter);
    filter_text="";
}

void System::set_filter(const std::vector<int> &ind)
{
    if(atoms.size()) throw PterosError("Filter could be set to empty system only!");
    filter = ind;
    sort_and_remove_duplicates(filter);
    filter_text="";
}

void System::set_filter(vector<int>::iterator it1, vector<int>::iterator it2)
{
    if(atoms.size()) throw PterosError("Filter could be set to empty system only!");
    copy(it1,it2,back_inserter(filter));
    sort_and_remove_duplicates(filter);
    filter_text="";
}

// Destructor of the system class
System::~System() {}

int System::frame_dup(int fr){
    if(fr<0 || fr>=traj.size())
        throw PterosError("Invalid frame {} for duplication!",fr);
    traj.push_back(traj[fr]);
    return traj.size()-1;
}

void System::frame_copy(int fr1, int fr2){
    if(fr1<0 || fr1>=traj.size() || fr2<0 || fr2>=traj.size())
        throw PterosError("Invalid frames {} and {} for copying!",fr1,fr2);
    traj[fr2] = traj[fr1];    
}

// Delete the range of frames. e = -1 is default
void System::frame_delete(int b, int e){
    int i;    

    if(e==-1) e = num_frames()-1;
    if(e<b || b<0 || e>num_frames()-1) throw PterosError("Invalid frame range {}:{} for deletion",b,e);

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
    if(traj.size()==0) LOG()->warn("All frames are deleted. All selections are now INVALID!");
}

void System::frame_swap(int fr1, int fr2)
{
    if(fr1<0 || fr1>=traj.size() || fr2<0 || fr2>=traj.size())
        throw PterosError("Invalid frames {} and {} for swapping!",fr1,fr2);
    std::swap(traj[fr1],traj[fr2]);
}

AtomHandler System::operator[](const std::pair<int, int> &ind_fr)
{
    return AtomHandler(*this,ind_fr.first,ind_fr.second);
}

/*
System::atom_iterator System::begin(int fr)
{
    return atom_iterator(*this,0,fr);
}

System::atom_iterator System::end(int fr)
{
    return atom_iterator(*this,num_atoms(),fr);
}
*/

void System::frame_append(const Frame& fr){
    traj.push_back(fr);
}

void System::assign_resindex(int start){
    if(start<0) start=0;

    int curres = atoms[start].resid;
    int curchain = atoms[start].chain;
    string curresname = atoms[start].resname;
    int cur = 0;
    if(start>0) cur = atoms[start].resindex;
    for(int i=start; i<atoms.size(); ++i){
        if( atoms[i].resid!=curres || atoms[i].chain!=curchain || atoms[i].resname!=curresname ){
            ++cur;
            curres = atoms[i].resid;
            curchain = atoms[i].chain;
            curresname = atoms[i].resname;
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
        if(atom(i).resindex == atom(j).resindex){
            return (i<j);
        } else {
            return (atom(i).resindex < atom(j).resindex);
        }
       }
    );
    // Now shuffle atoms and coordinates according to indexes
    vector<Atom> tmp(atoms); //temporary
    for(int i=0;i<ind.size();++i) atoms[i] = tmp[ind[i]];

    std::vector<Eigen::Vector3f> tmpv;
    for(int j=0; j<traj.size(); ++j){ // Over all frames
        tmpv = traj[j].coord; //temporary
        for(int i=0;i<ind.size();++i) traj[j].coord[i] = tmpv[ind[i]];

        if(traj[j].has_vel()){
            tmpv = traj[j].vel;
            for(int i=0;i<ind.size();++i) traj[j].vel[i] = tmpv[ind[i]];
        }

        if(traj[j].has_force()){
            tmpv = traj[j].force;
            for(int i=0;i<ind.size();++i) traj[j].force[i] = tmpv[ind[i]];
        }
    }
}

void System::clear_vel()
{
    for(int j=0; j<traj.size(); ++j) traj[j].vel.clear();
}

void System::clear_force()
{
    for(int j=0; j<traj.size(); ++j) traj[j].force.clear();
}

Selection System::atoms_dup(const vector<int>& ind){
    // Sanity check
    if(!ind.size()) throw PterosError("No atoms to duplicate!");
    for(int i=0; i<ind.size(); ++i){
        if(ind[i]<0 || ind[i]>atoms.size()-1)
            throw PterosError("Invalid index {} for atom duplication!", ind[i]);
    }

    // Duplicate atoms
    int first_added = atoms.size();
    int last_added = atoms.size()+ind.size()-1;
    // Prepare by increasing capacity of vectors
    atoms.reserve(atoms.size()+ind.size());
    for(int j=0; j<traj.size(); ++j){
        traj[j].coord.reserve(atoms.size()+ind.size());
        if(traj[j].has_vel()) traj[j].vel.reserve(atoms.size()+ind.size());
        if(traj[j].has_force()) traj[j].force.reserve(atoms.size()+ind.size());
    }

    // Now add atoms
    for(int i=0; i<ind.size(); ++i){
        // Add new atom
        atoms.push_back(atoms[ind[i]]);
        // Add new coordinate slot
        for(int j=0; j<traj.size(); ++j){
            traj[j].coord.push_back(traj[j].coord[ind[i]]);
            if(traj[j].has_vel()) traj[j].vel.push_back(traj[j].vel[ind[i]]);
            if(traj[j].has_force()) traj[j].force.push_back(traj[j].force[ind[i]]);
        }
    }


    return Selection(*this,first_added,last_added);
}

Selection System::atoms_add(const vector<Atom>& atm, const vector<Vector3f>& crd){
    // Sanity check
    if(!atm.size()) throw PterosError("No atoms to add!");
    if(atm.size()!=crd.size())
        throw PterosError("Number of coordinates {} doesn't mutch number of atoms {}",crd.size(),atm.size());

    int first_added = atoms.size();
    int last_added = atoms.size()+atm.size()-1;
    // Prepare by increasing capacity of vectors
    atoms.reserve(atoms.size()+atm.size());

    // If no frames add one
    if(traj.size()==0) traj.push_back(Frame());

    for(int j=0; j<traj.size(); ++j){
        traj[j].coord.reserve(atoms.size()+atm.size());
        if(traj[j].has_vel()) traj[j].vel.reserve(atoms.size()+atm.size());
        if(traj[j].has_force()) traj[j].force.reserve(atoms.size()+atm.size());
    }
    // Now add atoms
    for(int i=0; i<atm.size(); ++i){
        // Add new atom
        atoms.push_back(atm[i]);
        // Add new coordinate slot
        for(int j=0; j<traj.size(); ++j){
            traj[j].coord.push_back(crd[i]);
            if(traj[j].has_vel()) traj[j].coord.push_back(Vector3f::Zero());
            if(traj[j].has_force()) traj[j].coord.push_back(Vector3f::Zero());
        }
    }

    return Selection(*this,first_added,last_added);
}

void System::atoms_delete(const std::vector<int> &ind){
    int i,fr;

    // Sanity check
    if(!ind.size()) throw PterosError("No atoms to delete!");
    for(int i=0; i<ind.size(); ++i){
        if(ind[i]<0 || ind[i]>atoms.size()-1)
            throw PterosError("Index {} for atom is out of range (0:{})!",ind[i],atoms.size()-1);
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
    Frame tmpv;
    for(fr=0; fr<num_frames(); ++fr){
        // Make a copy
        tmpv = traj[fr];

        traj[fr].coord.clear();
        traj[fr].vel.clear();
        traj[fr].force.clear();
        for(i=0;i<tmp.size();++i){
            if(tmp[i].mass>=0){
                traj[fr].coord.push_back(tmpv.coord[i]);
                if(traj[fr].has_vel()) traj[fr].vel.push_back(tmpv.vel[i]);
                if(traj[fr].has_force()) traj[fr].force.push_back(tmpv.force[i]);
            }
        }
    }
}

void System::atom_move(int i, int j)
{
    // Sanity check
    if(i<0 || i>=num_atoms()) throw PterosError(format("Index of atom to move ({}}) is out of range ({}:{})!", i,0,num_atoms()));
    if(j<0 || j>=num_atoms()) throw PterosError(format("Target index to move ({}}) is out of range ({}:{}})!", j,0,num_atoms()));
    if(i==j) return; // Nothing to do

    // Move atom
    auto at = atom(i);

    if(i<j){
        for(int a=i+1; a<=j; ++a) atom(a-1) = atom(a);
        atom(j) = at;

        for(int fr=0; fr<num_frames(); ++fr){
            auto tmp = xyz(i,fr);
            for(int a=i+1; a<=j; ++a) xyz(a-1,fr) = xyz(a,fr);
            xyz(j,fr) = tmp;

            if(traj[fr].has_vel()){
                tmp = vel(i,fr);
                for(int a=i+1; a<=j; ++a) vel(a-1,fr) = vel(a,fr);
                vel(j,fr) = tmp;
            }

            if(traj[fr].has_force()){
                tmp = force(i,fr);
                for(int a=i+1; a<=j; ++a) force(a-1,fr) = force(a,fr);
                force(j,fr) = tmp;
            }
        }
    } else {
        for(int a=i-1; a>=j; --a) atom(a+1) = atom(a);
        atom(j) = at;

        for(int fr=0; fr<num_frames(); ++fr){
            auto tmp = xyz(i,fr);
            for(int a=i-1; a>=j; --a) xyz(a+1,fr) = xyz(a,fr);
            xyz(j,fr) = tmp;

            if(traj[fr].has_vel()){
                tmp = vel(i,fr);
                for(int a=i-1; a>=j; --a) vel(a+1,fr) = vel(a,fr);
                vel(j,fr) = tmp;
            }

            if(traj[fr].has_force()){
                tmp = force(i,fr);
                for(int a=i-1; a>=j; --a) force(a+1,fr) = force(a,fr);
                force(j,fr) = tmp;
            }
        }
    }
}

Selection System::atom_clone(int source)
{
    atoms_dup({source});
    atom_move(num_atoms()-1,source+1);
    return Selection(*this,source+1,source+1);
}

void System::atom_swap(int i, int j)
{
    if(i<0 || i>=num_atoms()) throw PterosError(format("Index of atom 1 to swap ({}}) is out of range ({}:{})!", i,0,num_atoms()));
    if(j<0 || j>=num_atoms()) throw PterosError(format("Index of atom 2 to swap ({}}) is out of range ({}:{}})!", j,0,num_atoms()));

    std::swap(atoms[i],atoms[j]);
    for(int fr=0; fr<num_frames(); ++fr){
        traj[fr].swap(i,j);
    }

}

Selection System::atom_add_1h(int target, int at1, int at2, int at3, float dist, bool pbc)
{
    Selection newat = atoms_dup({target});
    newat.name(0) = "H";
    newat.mass(0) = 1.0;
    newat.atomic_number(0) = 1;

    for(int fr=0; fr<num_frames(); ++fr){
        auto coor0 = xyz(target,fr);
        Vector3f coor1,coor2,coor3;

        if(pbc){
            coor1 = box(fr).closest_image(xyz(at1,fr),coor0);
            coor2 = box(fr).closest_image(xyz(at2,fr),coor0);
            coor3 = box(fr).closest_image(xyz(at3,fr),coor0);
        } else {
            coor1 = xyz(at1,fr);
            coor2 = xyz(at2,fr);
            coor3 = xyz(at3,fr);
        }

        // Coordinate of new atom
        Vector3f point = coor0-(coor1+coor2+coor3)/3.0;
        point = coor0-((point-coor1).cross(point-coor2)).normalized()*dist;

        newat.xyz(0,fr) = point;
    }

    return newat;
}

Selection System::atom_add_2h(int target, int at1, int at2, float dist, bool pbc)
{
    Selection newat1 = atoms_dup({target});
    newat1.name(0) = "H";
    newat1.mass(0) = 1.0;
    newat1.atomic_number(0) = 1;
    Selection newat2 = atoms_dup({newat1.index(0)});

    for(int fr=0; fr<num_frames(); ++fr){
        auto coor0 = xyz(target,fr);
        Vector3f coor1,coor2;

        if(pbc){
            coor1 = box(fr).closest_image(xyz(at1,fr),coor0);
            coor2 = box(fr).closest_image(xyz(at2,fr),coor0);
        } else {
            coor1 = xyz(at1,fr);
            coor2 = xyz(at2,fr);
        }

        // Coordinates of new atoms
        Vector3f up = (coor0-(coor1+coor2)/2.0).normalized()*dist;
        //Vector3f side = ((coor0-coor1).cross(coor0-coor2)).normalized();

        newat1.xyz(0,fr) = coor0+up;
        newat1.rotate(coor0,coor2-coor1,deg_to_rad(0.5*109.47));

        newat2.xyz(0,fr) = coor0+up;
        newat2.rotate(coor0,coor2-coor1,deg_to_rad(-0.5*109.47));
    }

    return Selection(*this,newat1.index(0),newat2.index(0));
}

Selection System::atom_add_3h(int target, int at1, float dist, bool pbc)
{
    Selection newat1 = atoms_dup({target});
    newat1.name(0) = "H";
    newat1.mass(0) = 1.0;
    newat1.atomic_number(0) = 1;
    Selection newat2 = atoms_dup({newat1.index(0),});
    Selection newat3 = atoms_dup({newat2.index(0),});

    for(int fr=0; fr<num_frames(); ++fr){
        auto coor0 = xyz(target,fr);
        Vector3f coor1;

        if(pbc){
            coor1 = box(fr).closest_image(xyz(at1,fr),coor0);
        } else {
            coor1 = xyz(at1,fr);
        }

        // Coordinates of new atoms
        Vector3f up = (coor0-coor1).normalized()/sqrt(24);
        Vector3f side = ((coor0-coor1).cross(Vector3f(1,0,0))).normalized()/sqrt(3);

        newat1.xyz(0,fr) = newat2.xyz(0,fr) = newat3.xyz(0,fr) = coor0+(up+side).normalized()*dist;

        auto m = rotation_transform(coor0, coor0-coor1, 2.0*M_PI/3.0);
        newat2.apply_transform(m);
        newat3.apply_transform(m*m);
    }
    return Selection(*this,newat1.index(0),newat3.index(0));
}

Selection System::append(const System &sys){
    //Sanity check
    if(num_frames()>0 && num_frames()!=sys.num_frames())
        throw PterosError("Can't merge systems with different number of frames ({} and {})!",num_frames(),sys.num_frames());
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
            traj[fr].time = sys.time(fr);
            traj[fr].box = sys.box(fr);
        }
        copy(sys.traj[fr].coord.begin(),sys.traj[fr].coord.end(),back_inserter(traj[fr].coord));
    }

    // Reassign resindex for added atoms
    assign_resindex(first_added-1);

    return Selection(*this,first_added,num_atoms()-1);
}

Selection System::append(const Selection &sel, bool current_frame){
    //Sanity check
    if(sel.size()==0) return Selection(*this); // In case of empty selection just exit
    if(!current_frame && num_frames()>0 && num_frames()!=sel.get_system()->num_frames())
        throw PterosError("Can't merge system with selection: different number of frames! ({} and {})",
                           num_frames(),sel.get_system()->num_frames());

    // If no frames create needed ammount
    bool transfer_time_box = false;    
    if(num_frames()==0){
        transfer_time_box = true;
        if(!current_frame){
            traj.resize(sel.get_system()->num_frames());
        } else {
            traj.resize(1);
        }
    }

    int first_added = num_atoms();    

    // Merge atoms
    atoms.reserve(atoms.size()+sel.size());
    for(int i=0;i<sel.size();++i) atoms.push_back(sel.atom(i));

    // Merge coordinates
    for(int fr=0; fr<num_frames(); ++fr){ // in system
        traj[fr].coord.reserve(atoms.size());
        if(transfer_time_box){
            if(!current_frame){
                // Transfer matching box and time
                traj[fr].time = sel.get_system()->time(fr);
                traj[fr].box = sel.get_system()->box(fr);
            } else {
                // Transfer box and time of current frame in sel
                traj[fr].time = sel.time();
                traj[fr].box = sel.box();
            }
        }
        for(int i=0;i<sel.size();++i){
            if(!current_frame){
                traj[fr].coord.push_back(sel.xyz(i,fr)); // Transfer matching frame coords
            } else {
                traj[fr].coord.push_back(sel.xyz(i)); // Transfer coords of current frame in sel
            }
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

Selection System::append(const AtomHandler &at)
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
        if(s.size()==0) LOG()->warn("Empty selection '{}' in rearrange will be ignored!", s.get_text());
        if(s.get_system()!=this) throw PterosError("Rearrange needs selections from the same system!");
    }
    // Overlap check
    if(check_selection_overlap(sel_vec)) throw PterosError("Selections for rearrange should not overlap!");

    System result;
    // Append all explicitly given selections to result
    for (auto &s: sel_vec) result.append(s);

    // Rest
    Selection rest(*this);
    rest.append(sel_vec);
    rest.invert();

    // append rest
    result.append(rest);

    *this = result;
}

void System::rearrange(const std::vector<string> &begin_strings, const std::vector<string> &end_strings)
{
    vector<Selection> begin_vec(begin_strings.size());
    for (int i=0; i<begin_strings.size(); ++i){
        begin_vec[i].modify(*this, begin_strings[i]);
    }

    vector<Selection> end_vec(end_strings.size());
    for (int i=0; i<end_strings.size(); ++i){
        end_vec[i].modify(*this, end_strings[i]);
    }

    rearrange(begin_vec, end_vec);
}

void System::rearrange(const std::vector<Selection> &begin_vec, const std::vector<Selection> &end_vec)
{
    // Index of middle
    int mid = begin_vec.size();

    // Combine vectors
    vector<Selection> vec;
    vec.insert(vec.end(),begin_vec.begin(),begin_vec.end());
    vec.insert(vec.end(),end_vec.begin(),end_vec.end());

    // Sanity check
    for(auto &s: vec){
        if(s.size()==0) LOG()->warn("Empty selection '{}' in rearrange will be ignored!", s.get_text());
        if(s.get_system()!=this) throw PterosError("Rearrange needs selections from the same system!");
    }

    // Overlap check
    if(check_selection_overlap(begin_vec)) throw PterosError("Selections for rearrange should not overlap!");

    // Compute rest
    Selection rest(*this);
    rest.append(vec);
    rest.invert();

    System result;
    // Append all starting selections to result
    for(int i=0;i<mid;++i) result.append(vec[i]);
    // Append rest
    result.append(rest);
    // Append all ending selections to result
    for(int i=mid;i<vec.size();++i) result.append(vec[i]);

    *this = result;
}

std::vector<std::pair<string, int> > System::rearrange_by_resname()
{
    // Find residue names
    set<string> rset;
    for(int i=0;i<num_atoms();++i) rset.insert(atoms[i].resname);

    // Sort residue names
    vector<string> rvec;
    std::copy(rset.begin(),rset.end(),back_inserter(rvec));
    std::sort(rvec.begin(),rvec.end());

    // Create selections and rearrange
    vector<Selection> selvec;
    for(auto& r: rvec) selvec.emplace_back(*this,"resname "+r);

    vector<pair<string,int>> res;
    for(int i=0;i<rvec.size();++i) res.emplace_back(rvec[i], selvec[i].num_residues());

    rearrange(selvec);

    return res;
}

void System::keep(const string &sel_str)
{
    Selection sel(*this,sel_str);
    keep(sel);
}

void System::keep(const Selection &sel)
{
    if(sel.get_system()!=this) throw PterosError("keep needs selection from the same system!");
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
    if(sel.get_system()!=this) throw PterosError("remove needs selection from the same system!");
    System tmp;
    tmp.append(~sel);
    *this = tmp;
    sel.clear();
}


void System::distribute(const Selection sel, Vector3i_const_ref ncopies, Matrix3f_const_ref shift)
{
    if(sel.get_system()!=this) throw PterosError("distribute needs selection from the same system!");

    Vector3f v;
    int last = num_atoms();
    int nc = 0;

    // How many atoms we are going to add
    int Nadd = (ncopies.prod()-1)*sel.size();

    // Allocate memory
    atoms.reserve(num_atoms()+Nadd);
    for(int fr=0; fr<num_frames(); ++fr){
        traj[fr].coord.reserve(atoms.size());
    }

    for(int x=0; x<ncopies(0); ++x){
        for(int y=0; y<ncopies(1); ++y){
            for(int z=0; z<ncopies(2); ++z){
                if(x>0 || y>0 || z>0){
                    // current shift
                    v = shift.col(0)*x + shift.col(1)*y + shift.col(2)*z;

                    for(int i=0;i<sel.size();++i){
                        atoms.push_back(sel.atom(i));
                    }
                    for(int fr=0; fr<num_frames(); ++fr){
                        for(int i=0;i<sel.size();++i) traj[fr].coord.push_back(sel.xyz(i,fr)+v);
                    }

                    ++nc;
                }
            }
        }
    }
    assign_resindex();
    //for(int i=0;i<num_atoms();++i) atoms[i].resid = atoms[i].resindex;
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

Selection System::select_all() const{
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
        v1 = box(fr).shortest_vector(xyz(i,fr),xyz(j,fr),pbc);
        v2 = box(fr).shortest_vector(xyz(k,fr),xyz(j,fr),pbc);
    } else {
        v1 = xyz(i,fr)-xyz(j,fr);
        v2 = xyz(k,fr)-xyz(j,fr);
    }
    return acos(v1.dot(v2)/(v1.norm()*v2.norm()));
}

float System::dihedral(int i, int j, int k, int l, int fr, Array3i_const_ref pbc) const
{
    Vector3f b1,b2,b3;
    if( (pbc!=0).any() ){
        Vector3f _i = xyz(i,fr);
        Vector3f _j = box(fr).closest_image(xyz(j,fr),_i,pbc);
        Vector3f _k = box(fr).closest_image(xyz(k,fr),_i,pbc);
        Vector3f _l = box(fr).closest_image(xyz(l,fr),_i,pbc);
        b1 = _j - _i;
        b2 = _k - _j;
        b3 = _l - _k;
    } else {
        b1 = xyz(j,fr)-xyz(i,fr);
        b2 = xyz(k,fr)-xyz(j,fr);
        b3 = xyz(l,fr)-xyz(k,fr);
    }    

    // Dihedral
    return atan2( ((b1.cross(b2)).cross(b2.cross(b3))).dot(b2/b2.norm()) ,
                  (b1.cross(b2)).dot(b2.cross(b3)) );
}

void System::wrap(int fr, Array3i_const_ref pbc){
    for(int i=0;i<num_atoms();++i){
        traj[fr].box.wrap_point(xyz(i,fr),pbc);
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


namespace pteros {

Vector2f get_energy_for_list(const vector<Vector2i>& pairs, const vector<float>& dist, const System& sys, vector<Vector2f>* pair_en){
    ForceField& ff = const_cast<System&>(sys).get_force_field();
    Vector2f e_total(0,0);

    if(pair_en) pair_en->resize(pairs.size());

    #pragma omp parallel
    {
        int at1,at2;
        Vector2f eloc(0,0);
        #pragma omp for nowait
        for(int i=0;i<pairs.size();++i){
            at1 = pairs[i](0);
            at2 = pairs[i](1);
            auto e = ff.pair_energy(at1, at2, dist[i],
                                sys.atom(at1).charge, sys.atom(at2).charge,
                                sys.atom(at1).type,   sys.atom(at2).type);
            if(pair_en) (*pair_en)[i] = e;
            eloc += e;
        }

        #pragma omp critical
        {
            e_total += eloc;
        }
    }
    return e_total;
}

//-------------------------------

void Frame::swap(int i, int j)
{
    std::swap(coord[i],coord[j]);
    if(has_vel()) std::swap(vel[i],vel[j]);
    if(has_force()) std::swap(force[i],force[j]);
}


}
