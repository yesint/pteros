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

//

#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/signals2.hpp>
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

    Energy_components(){
        total = 0.0;
        lj_14 = 0.0;
        q_14 = 0.0;
        lj_sr = 0.0;
        q_sr = 0.0;
    }
    /// Writes all energy components to string in the following order:
    /// total, lj_sr, lj_14, q_sr, q_14
    std::string to_str();
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
    float t;

    Frame(){        
        t = 0.0;
    }
};

//Forward declarations
class Selection;

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
    // Ensure correct 16-bytes-alignment for Eigen vectorization
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Constructors, operators and modification functions
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
    ~System();

    /// Append other system to this one
    void append(const System& sys);
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
    /// @name File IO
    /// @{

    /** Read structure, trajectory or topology from file.
    */
    // Skip functionality suggested by Raul Mera
    void load(std::string fname, int b=0, int e=-1, int skip = 0);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Operations with frames
    /// @{

    /// Duplicates given frame and adds it to the end of frame vector
    void frame_dup(int);

    /// Adds new frame to trajectory
    void frame_append(const Frame& fr);

    /// Copy coordinates from fr1 to fr2
    void frame_copy(int fr1, int fr2);

    /** Delete specified range of frames.
    *   If only @param b is supplied deletes all frames from b to the end.
    *   If only @param e is supplied deletes all frames from 0 to e
    */
    void frame_delete(int b = 0, int e = -1);

    /// Get read/write reference for given frame
    Frame& Frame_data(int fr);       
    /// @}
    ///

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Inline accessors
    /// @{

    /// Read/write access for periodic box for given frame
    inline Periodic_box& Box(int fr){
        return traj[fr].box;
    }

    /// Read/Write access to the time stamp of given frame
    inline float& Time(int fr){
    	return traj[fr].t;
    }

    /// Read/Write access for given coordinate of given frame
    inline Eigen::Vector3f& XYZ(int ind, int fr){
        return traj[fr].coord[ind];
    }
    /// @}


#ifndef NO_CPP11
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Secondary structure functions
    /// @{

    /// Determines secondary structure with DSSP algorithm and write detailed report to file
    void dssp(std::string fname) const;

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
#endif


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Manipulating atoms
    /// @{

    /// Adds new atoms, which are duplicates of existing ones by index
    void atoms_dup(const std::vector<int>& ind, Selection* res_sel = NULL);

    /// Adds new atoms from supplied vectors of atoms and coordinates
    void atoms_add(const std::vector<Atom>& atm,
                   const std::vector<Eigen::Vector3f>& crd,
                   Selection* res_sel = NULL);

    /// Delete the set of atoms
    void atoms_delete(const std::vector<int>& ind);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Periodicity-related functions
    /// @{

    /// Get distance between two atoms for given frame (periodic in given dimensions if needed).
    float distance(int i, int j, int fr, bool is_periodic = true, Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// Wrap all system to the periodic box for given frame
    void wrap_all(int fr, Vector3i_const_ref dims_to_wrap = Eigen::Vector3i::Ones());
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Energy functions
    /// @{

    /// Compute non-bond energy between two atoms
    /// The result is ADDED to e
    /// Intended mainly to be called from other functions, which take care of initializing e
    void add_non_bond_energy(Energy_components& e, int a1, int a2, int frame, bool is_periodic = true) const;

    /// Non-bond energy for given list of atom pairs
    Energy_components non_bond_energy(const std::vector<Eigen::Vector2i>& nlist, int fr) const;
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Utility functions
    /// @{

    /// Clears the system and prepares for loading completely new structure
    void clear(bool delete_selections = false);

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
            return NULL;
    }

    /// Assign unique resindexes
    /// This is usually done automatically upon loading a structure from file
    void assign_resindex();

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
