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


#pragma once

#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "pteros/core/atom.h"
#include "pteros/core/force_field.h"
#include "pteros/core/periodic_box.h"
#include "pteros/core/typedefs.h"
#include "pteros/core/atom_handler.h"


namespace pteros {

/// Definition of single trajectory frame.
/// Frames are stored in System class. They represent actual trajectory frames,
/// which are loaded from MD trajectories.
/// Coordinates are stored in nm as in Gromacs, not in Angstroms!
struct Frame {
    /// Coordinates of atoms
    std::vector<Eigen::Vector3f> coord;
    /// Velocities of atoms
    std::vector<Eigen::Vector3f> vel;
    /// Forces of atoms
    std::vector<Eigen::Vector3f> force;
    /// Periodic box
    PeriodicBox box;
    /// Timestamp
    float time;

    Frame(): time(0.0) {}

    bool has_vel() const { return !vel.empty(); }
    bool has_force() const { return !force.empty(); }

    /// Swap atoms i and j taking care of v and f
    void swap(int i, int j);
};

//====================================================================================

//Forward declarations
class Selection;
class FileHandler;
class FileContent;
class AtomHandler;

//====================================================================================

/**
*  The system of atoms.
*
*   The System is a container for atoms and their coordinates, which are typically
*   loaded from file. All properties of atoms, except the coordinates, are stored
*   in atoms vector. Coordinates are stored as a resizable vector of trajectory frames.
*   The system knows about selections associated with it and sends them
    signals if selections should adapt to the changes of coordinates of atom properties.
*   Copying and assignment of systems is allowed, but associated selections are
    not copied.
*/
class System {
    // System and Selection are friends because they are closely integrated.
    friend class Selection;    
    // Selection_parser must access internals of Selection
    friend class SelectionParser;
    // Needs an access for constructing the system in IO handlers
    friend class SystemBuilder;

public:    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Constructors and operators
    /// @{

    /// Default constructor
    System();

    /// Constructor creating system from file
    System(std::string fname);

    /// Copy constructor
    System(const System& other);

    /// Assignment operator
    System& operator=(const System& other);

    /// Constructing new System from Selection of other system
    System(const Selection& sel);

    /// Destructor
    virtual ~System();

    /// @}

    /// @name Adding, deleting and ordering groups of atoms
    /// These methods update resindexes automatically.
    /// @{

    /// Append other system to this one
    /// Returns selection corresponding to appended atoms
    Selection append(const System& sys);

    /// Append atoms from selection to this system.
    /// Selection may belong to any system, not necessary current.
    /// Returns selection corresponding to appended atoms.
    /// @param current_frame if true only the coordinates of current frame in sel are added to all frames of the system.
    Selection append(const Selection& sel, bool current_frame=false);

    /// Append single atom to this system
    /// Returns selection corresponding to appended atom
    Selection append(const Atom& at, Vector3f_const_ref coord);

    /** Append Atom_proxy object to this system
     Returns selection corresponding to appended atom.
     Usage:
     \code
     Selection sel(s,"name CA");
     System new_s;
     for(auto& at: sel){
        new_s.append(at);
     }
     \endcode
    */
    Selection append(const AtomHandler& at);

    /// Rearranges the atoms in the order of provided selection strings.
    /// Atom, which are not selected are appended at the end in their previous order.
    /// \note Selections should not overlap (exception is thrown if they are).
    void rearrange(const std::vector<std::string>& sel_strings);

    /// Rearranges the atoms in the order of provided selections.
    /// Atoms, which are not selected are appended as is at the end.
    /// \note Selections should not overlap (exception is thrown if they are).
    void rearrange(const std::vector<Selection>& sel_vec);

    /// Rearranges the atoms in the order of provided selection strings.
    /// @param begin_strings selection strings to begin with.
    /// @param end_strings   selection strings to end with.
    /// Atoms, which are not selected are placed as is between begin and end.
    /// \note Selections should not overlap (exception is thrown if they are).
    void rearrange(const std::vector<std::string>& begin_strings, const std::vector<std::string>& end_strings);

    /// Rearranges the atoms in the order of provided selections.
    /// @param begin_vec selections to begin with.
    /// @param end_vec   selections to end with.
    /// Atoms, which are not selected are placed as is between begin and end.
    /// \note Selections should not overlap (exception is thrown if they are).
    void rearrange(const std::vector<Selection>& begin_vec, const std::vector<Selection>& end_vec);

    /// Rearranges the atoms in the alphabetical order of residue names.
    /// Useful for liquid systems and membranes where each residue is a whole molecule.
    /// \warning If the molecules contain more than one resiude they will be broken!
    /// \returns vector of (residue name,number) pairs.
    std::vector<std::pair<std::string,int>> rearrange_by_resname();

    /// Keep only atoms given by selection string
    void keep(const std::string& sel_str);

    /// Keep only atoms from given selection
    void keep(const Selection& sel);

    /// Remove atoms given by selection string
    void remove(const std::string& sel_str);

    /// Remove atoms of given selection
    /// \warning Selection becomes invalid after that and is cleared!
    void remove(Selection& sel);

    /// Creates multiple copies of selection in the system and
    /// distributes them in a grid given by 3 translation vectors.
    /// Vectors are stored column-wise as shift.col(i)
    void distribute(const Selection sel, Vector3i_const_ref ncopies, Matrix3f_const_ref shift);

    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name General properties
    /// @{

    /// Returns the number of frames in selection
    inline int num_frames() const { return traj.size(); }

    /// Returns the number of atoms in selection
    inline int num_atoms() const { return atoms.size(); }
    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** @name Selecting atoms.
     These functions are just convenience adaptors for constructors of Selection class:
     \code
     // Conventional way:
     sys = System("structure.pdb");
     Selection sel(s,"name CA");

     // Shorter way:
     auto sel = sys.select("name CA");

     // Even shorter way using operator ():
     auto sel = sys("name CA");
     \endcode

     It also allows writing "one-liners" like this:
     \code
     System("file.pdb").select("name CA").write("ca.pdb");
     \endcode
    **/
    /// @{

    Selection select(std::string str, int fr = 0);
    Selection operator()(std::string str, int fr = 0);

    Selection select(int ind1, int ind2);
    Selection operator()(int ind1, int ind2);

    Selection select(const std::vector<int>& ind);
    Selection operator()(const std::vector<int>& ind);

    Selection select(std::vector<int>::iterator it1,
                     std::vector<int>::iterator it2);
    Selection operator()(std::vector<int>::iterator it1,
                         std::vector<int>::iterator it2);

    Selection select(const std::function<void(const System&,int,std::vector<int>&)>& callback, int fr = 0);
    Selection operator()(const std::function<void(const System&,int,std::vector<int>&)>& callback, int fr = 0);
    /// @}

    /// Convenience functions to select all
    /// @{
    Selection select_all() const;
    Selection operator()();
    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name File IO
    /// @{

    /** Read structure, trajectory or topology from file.
     @param b First frame to read
     @param e Last frame to read (-1 means up to the end of trajectory)
     @param skip Read only each skip frames
     @param on_frame Callback function, which takes pointer to the System as a
     first argument and the index of current stored frame as the second.
     If callback returns false loading stops.
    */
    void load(std::string fname,
              int b=0,
              int e=-1,
              int skip = 0,
              std::function<bool(System*,int)> on_frame = 0);

    /**
     * @brief Load data into System from the pre-opened file handler.
     * This is rather low-level method which provides
     * fine control over what should be read.
     * It can be called several times to read trajectory frames one by one
     * from the same pre-opened file.
     */
    bool load(const std::unique_ptr<FileHandler> &handler,
              FileContent what,
              std::function<bool(System*,int)> on_frame = 0);    

    void write(std::string fname, int b=-1,int e=-1) const;

    void write(const std::unique_ptr<FileHandler>& handler, FileContent what,int b=-1,int e=-1) const;


    /// Load Gromacs .ndx file and crease selections acording to it from existing system
    /// Returns a vector of pairs {name,Selection}
    std::vector<std::pair<std::string, Selection> > load_gromacs_ndx(std::string fname);

    /// @name Input filtering
    /** Filters narrow set of atoms and coordinates which are loaded from data files. Only atoms
        specified by filter are kept in the system.
        They follow the same rules as ordinary selections but can't be coordinate-dependent
        and can't be set by callbacs.
        \note Filters do not make loading faster. In fact it becomes slower
        (up to 2 times for very large frames) and consumes
        twise as much memory \e during loading. However when loading finishes the system will
        contain only selected atoms. If many frames are stored in the system overall
        memory consumption will be much smaller.

        Filter could only be set to empty system (it throws an error otherwise),
        thus the way of using them is the following:
        \code
        // Create empty system
        System sys;
        // Set filter
        sys.set_filter("name CA");
        // All calls to load will use filter now
        sys.load("some_protein.pdb");
        sys.load("trajectory.xtc");
        \endcode

        \warning Filters could not be used when loading topologies! load() throws an error if
        attempting to load topology with filter.
    */
    /// @{
    void set_filter(std::string str);
    void set_filter(int ind1, int ind2);
    void set_filter(const std::vector<int>& ind);
    void set_filter(std::vector<int>::iterator it1,
                    std::vector<int>::iterator it2);
    /// @}

    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Operations with frames
    /// @{

    /// Duplicates given frame and adds it to the end of frame vector
    int frame_dup(int fr);

    /// Appends provided frame to trajectory
    void frame_append(const Frame& fr);

    /// Copy all frame data from fr1 to fr2. Fr2 is overwritten!
    void frame_copy(int fr1, int fr2);

    /** Delete specified range of frames.
    *   If only @param b is supplied deletes all frames from b to the end.
    *   If only @param e is supplied deletes all frames from 0 to e
    */
    void frame_delete(int b = 0, int e = -1);

    /// Swaps two specified frames
    void frame_swap(int fr1, int fr2);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Inline accessors
    /// @{

    /// Read/write access for periodic box for given frame
    inline PeriodicBox& box(int fr=0){ return traj[fr].box; }

    /// Read only access for periodic box for given frame
    inline const PeriodicBox& box(int fr=0) const { return traj[fr].box; }

    /// Read/Write access to the time stamp of given frame
    inline float& time(int fr=0){ return traj[fr].time; }

    /// Read only access to the time stamp of given frame
    inline const float& time(int fr=0) const { return traj[fr].time; }

    /// Read/Write access for given coordinate of given frame
    inline Eigen::Vector3f& xyz(int ind, int fr=0){ return traj[fr].coord[ind]; }

    /// Read only access for given coordinate of given frame
    inline const Eigen::Vector3f& xyz(int ind, int fr=0) const { return traj[fr].coord[ind]; }

    /// Read/Write access for given velocity of given frame
    inline Eigen::Vector3f& vel(int ind, int fr=0){ return traj[fr].vel[ind]; }

    /// Read only access for given velocity of given frame
    inline const Eigen::Vector3f& vel(int ind, int fr=0) const { return traj[fr].vel[ind]; }

    /// Read/Write access for given force of given frame
    inline Eigen::Vector3f& force(int ind, int fr=0){ return traj[fr].force[ind]; }

    /// Read only access for given force of given frame
    inline const Eigen::Vector3f& force(int ind, int fr=0) const { return traj[fr].force[ind]; }

    /// Read/Write access for given atom
    inline Atom& atom(int ind) { return atoms[ind]; }

    /// Read only access for given atom
    inline const Atom& atom(int ind) const { return atoms[ind]; }

    /// Get read/write reference for given frame
    inline Frame& frame(int fr){ return traj[fr]; }

    /// Get read only reference for given frame
    inline const Frame& frame(int fr) const { return traj[fr]; }

    /// @}

    /// Iterators and indexing
    /// {@
    AtomIterator atoms_begin(){ return atoms.begin(); }
    AtomIterator atoms_end(){ return atoms.end(); }
    std::vector<Frame>::iterator traj_begin(){ return traj.begin(); }
    std::vector<Frame>::iterator traj_end(){ return traj.end(); }

    AtomHandler operator[](const std::pair<int,int>& ind_fr);

    /*
    class atom_iterator;
    /// Begin iterator. Takes frame number.
    atom_iterator begin(int fr);
    /// End iterator. Takes frame number.
    atom_iterator end(int fr);
    */

    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Manipulating sets of atoms by indexes.
    /// These methods <b>do not</b> update resindexes automatically.
    /// @{

    /// Adds new atoms, which are duplicates of existing ones by index. Atoms are placed at the end of the system.
    Selection atoms_dup(const std::vector<int>& ind);

    /// Adds new atoms from supplied vectors of atoms and coordinates. Atoms are placed at the end of the system.
    Selection atoms_add(const std::vector<Atom>& atm,
                   const std::vector<Eigen::Vector3f>& crd);

    /// Delete the set of atoms by indexes
    void atoms_delete(const std::vector<int>& ind);

    /// Move atom i to other position. Atom is inserted instead of atom j, shift is performed towards previous position of i.
    void atom_move(int i, int j);

    /// Duplicate single target atom and puts it immediately after the source.
    /// If using many times this procedure is slower than doing atoms_dup() and than assign_resindex().
    Selection atom_clone(int source);

    /// Swap two atoms
    void atom_swap(int i, int j);

    /// Add single hydrogen to target atom according to positions of 3 adjacent atoms assuming tetrahedral symmetry.
    Selection atom_add_1h(int target, int at1, int at2, int at3, float dist=0.109, bool pbc=true);

    /// Add two hydrogens to target atom according to positions of 2 adjacent atoms assuming tetrahedral symmetry.
    Selection atom_add_2h(int target, int at1, int at2, float dist=0.109, bool pbc=true);

    /// Add 3 hydrogens to target atom according to positions of single adjacent atom assuming tetrahedral symmetry.
    Selection atom_add_3h(int target, int at1, float dist=0.109, bool pbc=true);

    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Periodicity-related functions
    /// @{

    /// Wrap all system to the periodic box for given frame
    void wrap(int fr, Array3i_const_ref pbc = fullPBC);

    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Measuring functions
    /// @{

    /// Get distance between two atoms for given frame (periodic in given dimensions if needed).
    float distance(int i, int j, int fr, Array3i_const_ref pbc = fullPBC) const;

    /// Get angle in degrees between three atoms for given frame (periodic in given dimensions if needed).
    float angle(int i, int j, int k, int fr, Array3i_const_ref pbc = fullPBC) const;

    /// Get dihedral angle in degrees between three atoms for given frame (periodic in given dimensions if needed).
    float dihedral(int i, int j, int k, int l, int fr, Array3i_const_ref pbc = fullPBC) const;

    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Utility functions
    /// @{

    /// Clears the system and prepares for loading completely new structure
    void clear();

    bool force_field_ready(){return force_field.ready;}

    /// Returns internal Force_field object
    ForceField& get_force_field(){
        return force_field;
    }

    /// Assign unique resindexes
    /// This is usually done automatically upon loading a structure from file
    void assign_resindex(int start=0);

    /// Sorts atoms by resindex arranging atoms with the same resindexes
    /// into contigous pieces. Could be called after atom additions or duplications.
    void sort_by_resindex();

    /// Delete velocities from all frames.
    /// Does nothing if no velocities.
    void clear_vel();

    /// Delete forces from all frames.
    /// Does nothing if no forces.
    void clear_force();

    /// @}

protected:

    // Holds all atom attributes except the coordinates
    std::vector<Atom>  atoms;

    // Coordinates for any number of frames
    std::vector<Frame> traj;

    // Force field parameters
    ForceField force_field;

    // Indexes for filtering
    std::vector<int> filter;
    // Filter selection text for text-based filters
    std::string filter_text;

    void filter_atoms();
    void filter_coord(int fr);
};

//====================================================================================

/// Low level energy evaluation function
/// Returns total energy
/// Individual pair energies could be returned to pair_en if provided.
Eigen::Vector2f get_energy_for_list(const std::vector<Eigen::Vector2i>& pairs,
                                    const std::vector<float>& dist,
                                    const System& sys,
                                    std::vector<Eigen::Vector2f>* pair_en=nullptr);

//====================================================================================

/// Random-access forward iterator for Selection
/*
class System::atom_iterator {
public:
    using value_type = AtomHandler;
    using difference_type = size_t;
    using pointer = AtomHandler*;
    using reference = AtomHandler&;
    using iterator_category = std::random_access_iterator_tag;

    atom_iterator(const System& sys, int i, int fr): ind(i) {
        proxy.set(sys,ind,fr);
    }
    atom_iterator operator++(int junk) { atom_iterator tmp = *this; ++ind; return tmp; }
    atom_iterator& operator++() { ++ind; proxy.advance(); return *this; }
    atom_iterator& operator+(int i) {ind+=i; proxy.advance(i); return *this;}
    atom_iterator& operator-(int i) {ind-=i; proxy.advance(-i); return *this;}
    reference operator*() { return proxy; }
    pointer   operator->() { return &proxy; }
    bool operator==(const atom_iterator& rhs) { return ind == rhs.ind; }
    bool operator!=(const atom_iterator& rhs) { return !(*this==rhs); }

private:
    int ind;
    AtomHandler proxy;
};
*/

} // namespace


