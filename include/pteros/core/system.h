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

#ifndef SYSTEM_H
#define SYSTEM_H

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

namespace pteros {

/// Components of the non-bond energy
struct Energy_components {
    /// Total energy
    float total;
    /// Lenard-Jones energy of 1-4 pairs
    float lj_14;
    /// Coloumb energy of 1-4 pairs
    float q_14;
    /// Short-range Lenard-Jones energy (within cut-off)
    float lj_sr;
    /// Short-range Coloumb energy (within cut-off)
    float q_sr;

    Energy_components(): total(0.0), lj_14(0.0), q_14(0.0), lj_sr(0.0), q_sr(0.0) {}

    /// Writes all energy components to string in the following order:
    /// total, lj_sr, lj_14, q_sr, q_14
    std::string to_str();

    /// Addition of energies
    Energy_components operator+(const Energy_components& other);
    Energy_components& operator+=(const Energy_components& other);
};

/// Definition of single trajectory frame.
/// Frames are stored in System class. They represent actual trajectory frames,
/// which are loaded from MD trajectories.
/// Coordinates are stored in nm as in Gromacs, not in Angstroms!
struct Frame {
    /// Coordinates of atoms
    std::vector<Eigen::Vector3f> coord;
    /// Periodic box
    Periodic_box box;
    /// Timestamp
    float time;

    Frame(): time(0.0) {}
};

//Forward declarations
class Selection;
class Atom_proxy;
class Mol_file;
typedef unsigned short Mol_file_content;

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
    friend class Selection_parser;
    // Mol_file needs an access too
    friend class Mol_file;
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
    System& operator=(System other);

    /// Destructor
    virtual ~System();

    /// @}

    /// @name Adding, deleting and ordering groups of atoms
    /// These methods update resindexes automatically.
    /// @{

    /// Append other system to this one
    /// Returns selection corresponding to appended atoms
    Selection append(const System& sys);

    /// Append atoms from selection to this system
    /// Returns selection corresponding to appended atoms
    Selection append(const Selection& sel);

    /// Append single atom to this system
    /// Returns selection corresponding to appended atom
    Selection append(const Atom& at, const Vector3f_const_ref coord);

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
    Selection append(const Atom_proxy& at);

    /// Rearranges the atoms in the order of provided selection strings.
    /// Atom, which are not selected are appended at the end in their previous order.
    /// \note Selections should not overlap (exception is thrown if they are).
    void rearrange(const std::vector<std::string>& sel_strings);

    /// Rearranges the atoms in the order of provided selections.
    /// Atom, which are not selected are appended at the end in their previous order.
    /// \note Selections should not overlap (exception is thrown if they are).
    void rearrange(const std::vector<Selection>& sel_vec);

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
    /// distributes them in a grid
    void distribute(const Selection sel, Vector3i_const_ref ncopies, Vector3f_const_ref shift);

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

    Selection select(std::string str);
    Selection operator()(std::string str);

    Selection select(int ind1, int ind2);
    Selection operator()(int ind1, int ind2);

    Selection select(const std::vector<int>& ind);
    Selection operator()(const std::vector<int>& ind);

    Selection select(std::vector<int>::iterator it1,
                     std::vector<int>::iterator it2);
    Selection operator()(std::vector<int>::iterator it1,
                         std::vector<int>::iterator it2);

    Selection select(const std::function<void(const System&,int,std::vector<int>&)>& callback);
    Selection operator()(const std::function<void(const System&,int,std::vector<int>&)>& callback);

    /// Convenience function to select all
    Selection select_all();
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
    // Skip functionality suggested by Raul Mera
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
     * @param handler
     * @param what
     * @param b
     * @param e
     * @param skip
     * @param on_frame
     */
    void load(const std::unique_ptr<Mol_file> &handler,
              Mol_file_content what,
              int b=0,
              int e=-1,
              int skip = 0,
              std::function<bool(System*,int)> on_frame = 0);
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
    inline Periodic_box& Box(int fr){
        return traj[fr].box;
    }

    /// Read only access for periodic box for given frame
    inline const Periodic_box& Box(int fr) const {
        return traj[fr].box;
    }

    /// Read/Write access to the time stamp of given frame
    inline float& Time(int fr){
        return traj[fr].time;
    }

    /// Read only access to the time stamp of given frame
    inline const float& Time(int fr) const {
        return traj[fr].time;
    }

    /// Read/Write access for given coordinate of given frame
    inline Eigen::Vector3f& XYZ(int ind, int fr){
        return traj[fr].coord[ind];
    }

    /// Read only access for given coordinate of given frame
    inline const Eigen::Vector3f& XYZ(int ind, int fr) const {
        return traj[fr].coord[ind];
    }

    /// Read/Write access for given atom
    inline Atom& Atom_data(int ind) {
        return atoms[ind];
    }

    /// Read only access for given atom
    inline const Atom& Atom_data(int ind) const {
        return atoms[ind];
    }

    /// Get read/write reference for given frame
    inline Frame& Frame_data(int fr){
        return traj[fr];
    }

    /// Get read only reference for given frame
    inline const Frame& Frame_data(int fr) const {
        return traj[fr];
    }

    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Secondary structure functions
    /// @{

    /// Determines secondary structure with DSSP algorithm and writes detailed report to file
    void dssp(std::string fname) const;

    /// Determines secondary structure with DSSP algorithm and writes detailed report to stream
    void dssp(std::ostream& os) const;

    /**
     * @brief Determines secondary structure with DSSP algorithm and return it as a code string
     * @return Code string
     * The code is the same as in DSSP:
        alphahelix:	'H'
        betabridge:	'B'
        strand:		'E'
        helix_3:	'G'
        helix_5:	'I'
        turn:		'T'
        bend:		'S'
        loop:		' '
     */
    std::string dssp() const;
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Manipulating sets of atoms by indexes.
    /// These methods <b>do not</b> update resindexes automatically.
    /// @{

    /// Adds new atoms, which are duplicates of existing ones by index
    Selection atoms_dup(const std::vector<int>& ind);

    /// Adds new atoms from supplied vectors of atoms and coordinates
    Selection atoms_add(const std::vector<Atom>& atm,
                   const std::vector<Eigen::Vector3f>& crd);

    /// Delete the set of atoms by indexes
    void atoms_delete(const std::vector<int>& ind);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Periodicity-related functions
    /// @{

    /// Wrap all system to the periodic box for given frame
    void wrap_all(int fr, Vector3i_const_ref dims_to_wrap = Eigen::Vector3i::Ones());

    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Measuring functions
    /// @{

    /// Get distance between two atoms for given frame (periodic in given dimensions if needed).
    float distance(int i, int j, int fr, bool is_periodic = true,
                   Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// Get angle in degrees between three atoms for given frame (periodic in given dimensions if needed).
    float angle(int i, int j, int k, int fr, bool is_periodic = true,
                Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// Get dihedral angle in degrees between three atoms for given frame (periodic in given dimensions if needed).
    float dihedral(int i, int j, int k, int l, int fr, bool is_periodic = true,
                Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Energy functions
    /// @{

    /// Compute non-bond energy between two atoms        
    Energy_components non_bond_energy(int a1, int a2, int frame, bool is_periodic = true) const;

    /// Non-bond energy for given list of atom pairs
    Energy_components non_bond_energy(const std::vector<Eigen::Vector2i>& nlist, int fr, bool is_periodic=true) const;

    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Utility functions
    /// @{

    /// Clears the system and prepares for loading completely new structure
    void clear();

    /// Returns true if the force field is set up properly and is able to compute energies
    bool force_field_ready(){
        return force_field.ready;
    }

    /// Returns pointer to internal Force_field object for direct manipulation
    /// or NULL if force field is not ready
    Force_field* get_force_field(){
        if(force_field.ready)
            return &force_field;
        else
            return nullptr;
    }

    /// Assign unique resindexes
    /// This is usually done automatically upon loading a structure from file
    void assign_resindex(int start=0);

    /// Sorts atoms by resindex arranging atoms with the same resindexes
    /// into contigous pieces. Could be called after atom additions or duplications.
    void sort_by_resindex();

    /// @}

protected:

    /// Holds all atom attributes except the coordinates
    std::vector<Atom>  atoms;

    /// Coordinates for any number of frames
    std::vector<Frame> traj;

    /// Force field parameters
    Force_field force_field;

    /// Supplementary function to check if last added frame contains same number
    /// of atoms as topology
    void check_num_atoms_in_last_frame();
};

}
#endif /* SYSTEM_H */
