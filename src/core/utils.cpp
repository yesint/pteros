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





#include <iostream>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include "pteros/core/mol_file.h"
// DSSP
#include "pteros_dssp_wrapper.h"


using namespace std;
using namespace pteros;
using namespace Eigen;

// Base constructor of the system class
System::System() {

}

// Construnt system from file
System::System(string fname) {
    clear();
    load(fname);
}

System::System(const System& other){
    clear();
    atoms = other.atoms;
    traj = other.traj;
    force_field = other.force_field;
}

System& System::operator=(System other){
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
}

void System::check_num_atoms_in_last_frame(){
    if(Frame_data(num_frames()-1).coord.size()!=num_atoms())
        throw Pteros_error("File contains "
                           +to_string(Frame_data(num_frames()-1).coord.size())
                           +" atoms while the system has "
                           +to_string(num_atoms())
                           );
}

// Load structure or trajectory
void System::load(string fname, int b, int e, int skip, std::function<bool(System*,int)> on_frame){
    // Create an IO handler
    auto f = Mol_file::open(fname,'r');

    int num_stored = 0;    

    if(num_atoms()==0){
        // We don't have atoms yet, so we will read everything possible except trajectory
        auto c = f->get_content_type();

        if(c & MFC_COORD && !(c & MFC_TRAJ)){
            // We have single frame. Read it directly here
            Frame fr;
            frame_append(fr);
            f->read(this, &Frame_data(num_frames()-1), c);
            check_num_atoms_in_last_frame();
            ++num_stored;

            assign_resindex();

            // Call a callback if asked
            if(on_frame) on_frame(this,num_frames()-1);            
            // And now we should just exit
            return;
        } else {
            // We have not a single frame, so read only atoms here
            // Clear flags for trajectory and coordinates            
            c &= ~(MFC_COORD | MFC_TRAJ);
            f->read(this, nullptr, c);
            assign_resindex();
        }        
    }

    if(num_atoms()>0){
        // We have atoms already, so read only coordinates
        auto c = f->get_content_type();
        if(!(c & MFC_COORD) && !(c & MFC_TRAJ))
            throw Pteros_error("File reader for file '"+fname
                               +"' is not capable of appending frames to the system!");

        // Check if we can read trajectory
        if(c & MFC_TRAJ){
            // Sanity check for frame range
            if((e<b && e!=-1)|| b<0)
                throw Pteros_error("Invalid frame range for reading!");

            int cur = 0; // This holds real frame index in trajectory

            // Skip frames if needed
            if(b>0){
                cout << "Skipping " << b << " frames..." << endl;
                Frame skip_fr;
                for(int i=0;i<b;++i){
                    f->read(nullptr, &skip_fr, MFC_TRAJ);
                    cur++;
                }
            }

            int first = num_frames(); // Remember start

            cout << "Reading..."<<endl;

            int actually_read = 0;

            bool callback_ok = true;

            while(true){
                // End frame reached?
                if(cur==e && e!=-1) break;

                // Append new frame where the data will go
                Frame fr;
                frame_append(fr);
                // Try to read into it
                bool ok = f->read(nullptr, &Frame_data(num_frames()-1), MFC_TRAJ);
                if(!ok){
                    frame_delete(num_frames()-1); // Remove last frame - it's invalid
                    break; // end of trajectory
                }

                check_num_atoms_in_last_frame();

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
        } else if(c & MFC_COORD && !(c & MFC_TOP)) {
            // File contains single frame and no topology
            // Append new frame where the data will go
            Frame fr;
            frame_append(fr);
            // Read it
            f->read(nullptr, &Frame_data(num_frames()-1), MFC_COORD);
            check_num_atoms_in_last_frame();
            ++num_stored;
            // Call a callback if asked
            if(on_frame) on_frame(this,num_frames()-1);
        } else if(f->get_content_type() & MFC_TOP) {
            // This is topology file, read only topology                  
            f->read(this, nullptr, MFC_TOP);
        }
    }

    cout << "Accepted " << num_stored << " frames. Now " << num_frames() << " frames in the System" << endl;
}

// Destructor of the system class
System::~System() {}

int System::frame_dup(int fr){
    if(fr<0 || fr>=traj.size())
    	throw Pteros_error("Invalid frame for duplication!");
    traj.push_back(traj[fr]);
    return traj.size()-1;
}

void System::frame_copy(int fr1, int fr2){
    if(fr1<0 || fr1>=traj.size() || fr2<0 || fr2>=traj.size())
    	throw Pteros_error("Invalid frame for copying!");
    traj[fr2] = traj[fr1];    
}

// Delete the range of frames. e = -1 is default
void System::frame_delete(int b, int e){
    int i;    

    if(e==-1) e = num_frames()-1;
    if(e<b || b<0 || e>num_frames()-1) throw Pteros_error("Invalid frame range for deletion");    

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
    if(traj.size()==0) cout << "All frames are deleted. All selections are now INVALID!";
}

void System::frame_swap(int fr1, int fr2)
{
    if(fr1<0 || fr1>=traj.size() || fr2<0 || fr2>=traj.size())
        throw Pteros_error("Invalid frame for swapping!");
    std::swap(traj[fr1],traj[fr2]);
}

void System::frame_append(const Frame& fr){
    traj.push_back(fr);
}

void System::assign_resindex(int start){
    if(start<0) start=0;
    //cout << "Assigning resindex..." << endl;
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
    sort(ind.begin(),ind.end(), [this](int i, int j){
                                    return (Atom_data(i).resindex < Atom_data(j).resindex);
                                } );
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
            throw Pteros_error("Invalid index for atom duplication!");
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
    if(atm.size()!=crd.size()) throw Pteros_error("Wrong number of coordinates for adding atoms!");

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
            throw Pteros_error("Invalid index for atom deletion!");
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

Selection System::append(const System &sys){
    //Sanity check
    if(num_frames()>0 && num_frames()!=sys.num_frames()) throw Pteros_error("Can't merge systems with different number of frames!");
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
    if(num_frames()>0 && num_frames()!=sel.get_system()->num_frames()) throw Pteros_error("Can't merge systems with different number of frames!");
    // If no frames create needed ammount
    bool transfer_time_box = false;
    if(num_frames()==0){
        transfer_time_box = true;
        traj.resize(sel.get_system()->num_frames());        
    }

    int first_added = num_atoms();    

    // Merge atoms
    //atoms.reserve(atoms.size()+sel.size());
    for(int i=0;i<sel.size();++i) atoms.push_back(sel.Atom_data(i));
    // Merge coordinates    
    for(int fr=0; fr<num_frames(); ++fr){
        //traj[fr].coord.reserve(atoms.size()+sel.size());
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

Selection System::append(const Atom &at, const Vector3f_const_ref coord)
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
    append(at.Atom_data(),at.XYZ());
}

void System::rearrange(const std::vector<string> &sel_strings){
    vector<Selection> sel_vec(sel_strings.size());
    for (int i=0; i<sel_strings.size(); ++i){
        sel_vec[i].modify(*this, sel_strings[i]);
    }

    rearrange(sel_vec);
}

void System::rearrange(const std::vector<Selection> &sel_vec){
    // Empty selections check
    for(auto &s: sel_vec){
        if(s.size()==0) throw Pteros_error("Empty selections are not permitted in rearrange!");
    }
    // Overlap check
    vector<int> inters;
    for (int i=0; i<sel_vec.size()-1; ++i){
        for (int j=i+1; j<sel_vec.size(); ++j){
            set_intersection(sel_vec[i].index_begin(), sel_vec[i].index_end(),
                             sel_vec[j].index_begin(), sel_vec[j].index_end(),
                             back_inserter(inters));
            if (!inters.empty()) throw Pteros_error("Selections for rearrange should not overlap!");
        }
    }

    Selection rest(*this);
    for (auto &s: sel_vec) rest.append(s);
    rest = ~rest;

    System result;
    for (auto &s: sel_vec) result.append(s);
    result.append(rest);

    *this = result;
}

void System::keep(const string &sel_str)
{
    Selection sel(*this,sel_str);
    keep(sel);
}

void System::keep(const Selection &sel)
{
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
    System tmp;
    tmp.append(~sel);
    *this = tmp;
    sel.clear();
}

void System::distribute(const Selection sel, Vector3i_const_ref ncopies, Vector3f_const_ref shift)
{
    Selection tmp;
    for(int x=0; x<ncopies(0); ++x)
        for(int y=0; y<ncopies(1); ++y)
            for(int z=0; z<ncopies(2); ++z){
                if(x>0 || y>0 || z>0){
                    tmp = append(sel);                    
                    tmp.translate(Vector3f(shift(0)*x,shift(1)*y,shift(2)*z));
                }
            }
}

Selection System::select(string str){
    return Selection(*this,str);
}

Selection System::select(int ind1, int ind2){
    return Selection(*this,ind1,ind2);
}

Selection System::select(const std::vector<int> &ind){
    return Selection(*this,ind);
}

Selection System::select(std::vector<int>::iterator it1, std::vector<int>::iterator it2){
    return Selection(*this,it1,it2);
}

Selection System::select_all(){
    return Selection(*this,0,num_atoms()-1);
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

float System::distance(int i, int j, int fr, bool is_periodic, Vector3i_const_ref dims) const {
    if(is_periodic){
        return traj[fr].box.distance(traj[fr].coord[i], traj[fr].coord[j], dims);
    } else {
        return (traj[fr].coord[i] - traj[fr].coord[j]).norm();
    }
}

#define RAD_TO_DEG 57.295779513082320876798154814105

float System::angle(int i, int j, int k, int fr, bool is_periodic, Vector3i_const_ref dims) const
{
    Vector3f v1,v2;
    if(is_periodic){
        v1 = Box(fr).shortest_vector(XYZ(i,fr),XYZ(j,fr),dims);
        v2 = Box(fr).shortest_vector(XYZ(k,fr),XYZ(j,fr),dims);
    } else {
        v1 = XYZ(i,fr)-XYZ(j,fr);
        v2 = XYZ(k,fr)-XYZ(j,fr);
    }
    return acos(v1.dot(v2)/(v1.norm()*v2.norm())) * RAD_TO_DEG;
}

float System::dihedral(int i, int j, int k, int l, int fr, bool is_periodic, Vector3i_const_ref dims) const
{
    Vector3f b1,b2,b3;
    if(is_periodic){
        b1 = Box(fr).get_closest_image(XYZ(j,fr),XYZ(i,fr)) - XYZ(i,fr);
        b2 = Box(fr).get_closest_image(XYZ(k,fr),XYZ(i,fr)) -
               Box(fr).get_closest_image(XYZ(j,fr),XYZ(i,fr));
        b3 = Box(fr).get_closest_image(XYZ(l,fr),XYZ(i,fr)) -
               Box(fr).get_closest_image(XYZ(k,fr),XYZ(i,fr));
    } else {
        b1 = XYZ(j,fr)-XYZ(i,fr);
        b2 = XYZ(k,fr)-XYZ(j,fr);
        b3 = XYZ(l,fr)-XYZ(k,fr);
    }    

    // Dihedral
    return atan2( ((b1.cross(b2)).cross(b2.cross(b3))).dot(b2/b2.norm()) ,
                  (b1.cross(b2)).dot(b2.cross(b3)) ) * RAD_TO_DEG;
}

void System::wrap_all(int fr, Vector3i_const_ref dims_to_wrap){
    for(int i=0;i<num_atoms();++i){
        traj[fr].box.wrap_point(XYZ(i,fr),dims_to_wrap);
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

string Energy_components::to_str(){
    return    to_string(total) + " "
            + to_string(lj_sr) + " "
            + to_string(lj_14) + " "
            + to_string(q_sr) + " "
            + to_string(q_14);
}

Energy_components Energy_components::operator+(const Energy_components& other){
    Energy_components ret;
    ret.total = total+other.total;
    ret.lj_14 = lj_14+other.lj_14;
    ret.q_14 = q_14+other.q_14;
    ret.lj_sr = lj_sr+other.lj_sr;
    ret.q_sr = q_sr+other.q_sr;
    return ret;
}

Energy_components &Energy_components::operator+=(const Energy_components &other)
{
    total += other.total;
    lj_14 += other.lj_14;
    q_14 += other.q_14;
    lj_sr += other.lj_sr;
    q_sr += other.q_sr;
    return *this;
}

Energy_components System::non_bond_energy(int a1, int a2, int frame, bool is_periodic) const
{
    Energy_components e;

    // First check if this pair is not in exclusions
    if( force_field.exclusions[a1].count(a2) == 0 ){
        // Required at1 < at2
        int at1,at2;
        if(a1<a2){
            at1 = a1;
            at2 = a2;
        } else {
            at1 = a2;
            at2 = a1;
        }

        float e1,e2;

        int N = force_field.LJ14_interactions.size();

        float r = distance(at1,at2,frame,is_periodic);

        // Check if this is 1-4 pair
        std::unordered_map<int,int>::iterator it = const_cast<System&>(*this).force_field.LJ14_pairs.find(at1*N+at2);
        if( it == force_field.LJ14_pairs.end() ){
            // Normal, not 1-4
            e1 = force_field.LJ_kernel_ptr(force_field.LJ_C6(atoms[at1].type,atoms[at2].type),
                                   force_field.LJ_C12(atoms[at1].type,atoms[at2].type),
                                   r);
            e2 = force_field.coulomb_kernel_ptr(atoms[at1].charge,
                                       atoms[at2].charge,
                                       r);
            e.lj_sr += e1;
            e.q_sr += e2;
            e.total += (e1 + e2);
        } else {
            // 1-4
            e1 = force_field.LJ_kernel_ptr(force_field.LJ14_interactions[it->second](0),
                                   force_field.LJ14_interactions[it->second](1),
                                   r);
            e2 = force_field.coulomb_kernel_ptr(atoms[at1].charge,
                                       atoms[at2].charge,
                                       r)
                    * force_field.fudgeQQ;

            e.lj_14 = e1;
            e.q_14 = e2;
            e.total += (e1 + e2);
        }
    }

    return e;
}

Energy_components System::non_bond_energy(const std::vector<Eigen::Vector2i> &nlist, int fr, bool is_periodic) const
{
    Energy_components e;
    int n = nlist.size();

    #pragma omp parallel
    {
        Energy_components _e;
        #pragma omp for nowait
        for(int i=0;i<n;++i){
            _e += non_bond_energy(nlist[i](0),nlist[i](1),fr,is_periodic);
        }
        #pragma omp critical
        {
            e += _e;
        }
    }

    return e;
}

void System::dssp(string fname) const {
    ofstream f(fname.c_str());
    Selection sel(const_cast<System&>(*this),"all");
    dssp_wrapper(sel,f);
    f.close();
}

void System::dssp(ostream& os) const {
    Selection sel(const_cast<System&>(*this),"all");
    dssp_wrapper(sel,os);
}


string System::dssp() const{
    Selection sel(const_cast<System&>(*this),"all");
    return dssp_string(sel);
}




