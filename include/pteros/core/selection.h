/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2014, Semen Yesylevskyy
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

#ifndef SELECTION_H
#define SELECTION_H

#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/signals2.hpp>
#include "pteros/core/selection_parser.h"
#include "pteros/core/system.h"
#include <boost/shared_ptr.hpp>

namespace pteros {

// Forward declaration of frient class
class System;
class Selection_parser;

/** @brief Selection class.
*
*   Selections are key objects in Pteros. Technically speaking the selection
*   is just an array, which contains indexes of selected atoms in particlar system.
*   Selection does not hold the copies of the atoms or their coordinates, it
*   just points to them serving like a handy alias for certain subset of atoms.
*   Thus selections may overlap arbitrarily.
*   Selections are used to perform various operations on the group of selected atoms.
*   The changes become immediately visible to all other selections, which point to
*   some of changed atoms.
*   Each selection is bound to particular System. There are neither 'parentless' selection nor the
*   selections, which combine the atoms from different systems.
*   Selections are created using the syntax, which is very similar to those used in VMD.
*/
class Selection {
  /// System and Selection are friends because they are closely integrated.
  friend class System;
  friend class Selection_parser;

  public:
    /// Ensure correct 16-bytes-alignment for Eigen sse2 optimizations
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /** Main constructor.
    *   @param sys System pointed by this selection
    *   @param str Selection string
    *   if with_signal is true System with send signals to selection automatically
    */
    Selection(System& sys, std::string str);

    /** Constructor with delayed parsing.
    *   Associates selection with the system @param sys,
        but does not parse selection.
    *   Selection text should be passed later by overloaded << operator or by
    *   calling modify() function.
    *   if with_signal is true System with send signals to selection automatically
    */
    Selection(System& sys);

    /// Default constructor for absolutely empty selection.
    Selection();

    /** Constructor, which creates selection from the interval of indexes
        instead of selection string.
        It is much faster then parsing corresponding string, but is limited
        to contigous interval of indexes.
        if with_signal is true System with send signals to selection automatically
     */
    Selection(System& sys, int ind1, int ind2);

    /// Assignment operator
    Selection& operator=(Selection);    

    /// Copy constructor
    Selection(const Selection&);

    /// Equality operators
    bool operator==(const Selection &other) const;

    bool operator!=(const Selection &other) const {
        return !(*this == other);
    }

    /// Destructor
    ~Selection();

    /// Append another selection to this one. Acts like logical "OR" between selections.
    void append(Selection& sel);

    /// Append absolute index to selection
    void append(int ind);

    /// Modifies both system and string in selection.
    /// If with_signal is true System with send signals to selection automatically
    void modify(System& sys, std::string str);

    /** Modifies selection string in existing selection.
    *   @param str New value of selection text. Selection is re-parsed immediately with
    *   this new value.
    */
    void modify(std::string str);

    void modify(System& sys);

    /// Modifies both system and selection using the range of indexes.
    /// If with_signal is true System with send signals to selection automatically.
    void modify(System& sys, int ind1, int ind2);

    /// Modifies selection using the range of indexes
    void modify(int ind1, int ind2);

    /// Modifies selection using vector of indexes
    void modify(std::vector<int>& ind);

    /// Modifies selection using pair of iterators to index array
    void modify(std::vector<int>::iterator it1, std::vector<int>::iterator it2);

    /** Recomputes selection without re-parsing selection text.
    *   Only makes sense for coordinate-dependent selections when the coordinates change.
    *   Called authomatically for coordinate-dependent selections by set_frame()
    *   If selection is not coordinate-dependent does nothing.
    */
    void apply();

    /** Recomputes selection completely.
    *   May be used when new file is loaded into the system, or when atoms are
    *   created/deleted. Forces re-parsing of selection text.
    */
    void update();

    /// @name Frame functions
    /// @{

    /// Get current frame selection is pointing to
    int get_frame() const {return frame;}
    /// Set current frame for selection
    void set_frame(int fr);
    /// @}

    /// Get the size of selection
    int size() const {return index.size();}

    /** Clears selection and frees memory, but do not delete it.
    *   Selection remains registered in parent system, but is cleared from any data.
    *   Subsequent call of modify() may populate it again.
    */
    void clear();

    /// Selects each residue, which is references by selection.
    /// All selections for residues are placed into supplied vector.
    /// Handles multiple chains correctly
    /// If with_signal is true selections are created with automatic signalling from parent system
    void each_residue(std::vector<Selection>& sel) const;

    /// @name Get and set functions
    /// @{

    /// Get pointer to the system, which owns this selection
    System* get_system() const { return system; }
    /// Get selection text
    std::string get_text() const;
    /// Get vector of all indexes in selection
    std::vector<int> get_index() const { return index; }
    /// Get const iterator for begin of index
    std::vector<int>::const_iterator index_begin() const { return index.begin(); }
    /// Get const iterator for the end of index
    std::vector<int>::const_iterator index_end() const { return index.end(); }
    /// Get vector of all chains in selection
    std::vector<char> get_chain() const;
    /// Set chains from supplied vector.
    /// Its size must be the save as the size of selection.
    void set_chain(const std::vector<char>& data);
    /// Sets chain of all selected atoms to the same given value.
    void set_chain(char data);
    /// Get vector of unique chains in selection
    std::vector<char> get_unique_chain() const;
    /// Get vector of all resid's in selection
    /// This works correctly inside one chain only!
    /// For multiple chains resid's will overlap.
    std::vector<int> get_resid() const;
    /// Get vector of unique resid's in selection
    std::vector<int> get_unique_resid() const;
    /// Set resid's in selection from supplied vector.
    /// Its size must be the save as the size of selection.
    void set_resid(const std::vector<int>& data);
    /// Sets resid of all selected atoms to the same given value.
    void set_resid(int data);
    /// Get vector of all resindexes in selection. Resindexes are unique
    /// regardless the number of the chains.
    std::vector<int> get_resindex() const;
    /// Get vector of unique resindexes's in selection
    std::vector<int> get_unique_resindex() const;
    /// Get vector of all atom names in selection
    std::vector<std::string> get_name() const;
    /// Set atom names in selection from supplied vector.
    /// Its size must be the save as the size of selection.
    void set_name(const std::vector<std::string>& data);
    /// Sets atom names of all selected atoms to the same given value.
    void set_name(std::string& data);

    /// Get vector of all resnames in selection
    std::vector<std::string> get_resname() const;
    /// Set resnames in selection from supplied vector.
    /// Its size must be the save as the size of selection.
    void set_resname(const std::vector<std::string>& data);
    /// Sets resnames of all selected atoms to the same given value.
    void set_resname(std::string& data);

    /// Get coordinates of all atoms in this selection for current frame
    Eigen::MatrixXf get_xyz() const;
    void get_xyz(Eigen::MatrixXf& res) const;
    /// Set coordinates of this selection for current frame
    void set_xyz(const Eigen::MatrixXf& coord);
    /// Computes average structure over the range of frames
    Eigen::MatrixXf get_average(int b=0, int e=-1) const;
    /// Get masses of all atoms in selection
    std::vector<float> get_mass() const;
    /// Set atom masses in selection to the values from supplied vector.
    /// Its size must be the save as the size of selection.
    void set_mass(const std::vector<float> m);
    /// Sets masses of all selected atoms to the same given value.
    void set_mass(float data);
    /** Extracts X,Y,Z for given atom index for specified range of frames
        (gets trajectory of given atom).
    *   Result is returned as MatrixXf, where i-th column is an XYZ vector for frame i.    
    */
    Eigen::MatrixXf get_traj(int ind, int b=0, int e=-1) const;

    /// Get beta
    std::vector<float> get_beta() const;
    /// Set beta in selection to the values from supplied vector.
    /// Its size must be the save as the size of selection.
    void set_beta(std::vector<float>& data);
    /// Sets beta of all selected atoms to the same given value.
    void set_beta(float data);    
    /// @}

    /// @name Inquery functions
    /// @{

    /** Get the center of selection.
    * @param mass_weighted Use mass-weighting
    * @param periodic Account for periodic boundary conditions.
    * Please note that if the size of selection is larger than 1/2 of the box size in
    * any dimension you will get incorrect results if periodic is set to true.
    */
    Eigen::Vector3f center(bool mass_weighted = false, bool periodic = false) const;
    /// Get minimal and maximal coordinates in selection
    void minmax(Eigen::Vector3f& min, Eigen::Vector3f& max) const;

#ifdef USE_POWERSASA
    /// Get the SASA. Easy way - only returns SASA area of selection
    float sasa(float probe_r = 0.14);

    /// Get the SASA. Detailed way - returns area and computes volume and per-atom values
    float sasa(float probe_r = 0.14, float* total_volume = NULL,
               std::vector<float>* area_per_atom = NULL,
               std::vector<float>* volume_per_atom = NULL);
#endif
    /// @}

    /// @name Geometry transformation functions
    /// @{

    /// Translate selection by given vector
    void translate(const Eigen::Vector3f&);
    /// Rotate around given axis relative to cm
    /// @param axis Axis of rotation 0=X, 1=Y, 2=Z
    /// @param angle Rotation angle in radians
    void rotate(int axis, float angle);
    /// Rotate around given axis relative to given pivot
    /// @param axis Axis of rotation 0=X, 1=Y, 2=Z
    /// @param angle Rotation angle in radians
    /// @param pivot Rotation around this pivot
    void rotate(int axis, float angle, const Eigen::Vector3f& pivot);
    /// Rotate around given vector relative to given pivot
    /// @param direction Rotate around this vector
    /// @param angle Rotation angle in radians
    /// @param pivot Rotation around this pivot
    void rotate(const Eigen::Vector3f& direction, float angle, const Eigen::Vector3f& pivot);
    /// Rotation with the given 3x3 rotation matrix around point (0,0,0)
    void rotate(const Eigen::Matrix3f& m);
    /// Rotation by given angles around X, Y and Z with given pivot
    void rotate(const Eigen::Vector3f& angles, const Eigen::Vector3f& pivot);

    /// Wraps whole selection to the periodic box
    void wrap(const Eigen::Vector3i& dims = Eigen::Vector3i::Ones());

    /** Unwraps selection to make it whole if possible (without jumps over periodic box boundary).
     * The periodic center of mass is used as an anchor point.
     * Please note that if the size of selection is larger than 1/2 of the box size in
     * any dimension unwrap() will not work as expected and will not make selection "compact"!
    */
    void unwrap(const Eigen::Vector3i& dims = Eigen::Vector3i::Ones());

    /** Unwraps selection to make it whole (without jumps over periodic box boundary).
     * based on preserving all bonds. The maximal bond length is given by d.
     * This method works reliably in any case, but is much slower than unwrap()
     */
    void unwrap_bonds(float d = 0.2, const Eigen::Vector3i& dims = Eigen::Vector3i::Ones());

    /** Get transform for orienting selection by principal axes.
     * Please note that if the size of selection is larger than 1/2 of the box size in
     * any dimension you will get funny results if is_periodic is set to true.
     */
    Eigen::Affine3f principal_transform(bool is_periodic = false);

    /** Orient molecule by its principal axes.
     * The same as
     *
     * Eigen::Affine3f tr = sel.principal_transform();
     * sel.apply_transform(tr);
     */
    void principal_orient(bool is_periodic = false);
    /// @}

    /// @name Fitting and RMSD functions
    /// @{

    /// RMSD between current and another frame
    float rmsd(int fr);

    /// RMSD between two frames
    float rmsd(int fr1, int fr2);

    /// RMSD between two selections of the same size
    friend float rmsd(Selection& sel1, Selection& sel2);

    /** RMSD between two selections of the same size (for given frames)
    *   @param sel1 First selection
    *   @param fr1 Frame for first selection
    *   @param sel2 Second selection
    *   @param fr2 Frame for second selection
    */
    friend float rmsd(Selection& sel1, int fr1, Selection& sel2, int fr2);

    /// Fit two selection of the same size
    friend void fit(Selection& sel1, Selection& sel2);

    /// Fit all frames in the trajectory to reference frame
    void fit_trajectory(int ref_frame=0, int b=0, int e=-1);

    /// Returns fitting transformation for two given selections of the same size
    friend Eigen::Affine3f fit_transform(Selection&, Selection&);    

    /// Returns fit transformation between frames fr1 and fr2
    Eigen::Affine3f fit_transform(int fr1, int fr2);

    /// Fits frame fr1 to fr2
    void fit(int fr1, int fr2);

    /// Apply fitting transformation
    void apply_transform(Eigen::Affine3f& t);
    /// @}


    /// @name Energy functions
    /// @{

    /// Self-energy of selection
    Energy_components non_bond_energy() const;
    /// Non-bond energy between two selections
    friend Energy_components non_bond_energy(Selection& sel1, Selection& sel2, int fr);    

    /// @}

    /// @name File IO
    /// @{

    /** Write structure or trajectory for selection. Type is determined by extension.
    *   Frames from b to e are written.
    *   If b and e are not set they default to current frame
    */
    void write(std::string fname,int b=-1,int e=-1);
    /// @}


    /// @name Building functions
    /// @{
    /// Duplicate current selection in the parent system
    /// Resulting selection is born without signalling because it is intended for building purposes only
    void atoms_dup(Selection* res_sel = NULL);
    /// Delete all atoms of current selection from the parent system
    void atoms_delete();
    /// Creates multiple copies of selection in the parent system and
    /// distributes them in a grid
    void distribute(Eigen::Vector3i& ncopies, Eigen::Vector3f& shift);
    /// @}

    /// @name Signals handling
    /// @{
    /// Enable automatic signalling from system
    void enable_signals();
    /// Disable automatic signalling from system
    void disable_signals();
    /// Returns signalling state of selection
    bool signals_enabled() const;
    /// @}

    /// Split current selection into several selections according to
    /// the interatomic distances. Each resulting selection is a group
    /// of atoms connected by distances less than d
    void split_by_connectivity(float d, std::vector<Selection>& res);

    void split_by_residue(std::vector<Selection>& res);

    /// Computes the central momens of inertia and principal axes of inertia
    void inertia(Eigen::Vector3f& moments, Eigen::Matrix3f& axes, bool periodic = false) const;

    /// Computes radius of gyration for selection
    float gyration(bool periodic = false) const;

    /** @name Inlined utility functions.
    *   Used to access the properties of
    *   particular atom in selection. The ind passed to these functions is
    *   is the selection index, not the global system %index. I.e. passing 10 will
    *   extract the property of %atom with global %index index[10]
    *   <br>All these function except Index could be used as lvalue, which means that
    *   it is possible to assign value to them. For example X(10) = 3.14 will assign the
    *   value 3.14 to the atom with global index index[10].
    */
    /// @{

    /// Extracts X for current frame
    inline float& X(int ind){
        return system->traj[frame].coord[index[ind]](0);
    }

    /// Extracts X for given frame frame fr
    inline float& X(int ind, int fr){
        return system->traj[fr].coord[index[ind]](0);
    }

    /// Extracts Y for current frame
    inline float& Y(int ind){
        return system->traj[frame].coord[index[ind]](1);
    }

    /// Extracts Y for given frame frame fr
    inline float& Y(int ind, int fr){
    	return system->traj[fr].coord[index[ind]](1);
    }

    /// Extracts Z for current frame
    inline float& Z(int ind){
    	return system->traj[frame].coord[index[ind]](2);
    }

    /// Extracts Z for given frame frame fr
    inline float& Z(int ind, int fr){
    	return system->traj[fr].coord[index[ind]](2);
    }

    /// Extracts X,Y and Z for current frame
    inline Eigen::Vector3f& XYZ(int ind){
    	return system->traj[frame].coord[index[ind]];
    }

    /// Extracts X,Y and Z for given frame frame fr
    inline Eigen::Vector3f& XYZ(int ind, int fr){
    	return system->traj[fr].coord[index[ind]];
    }

    /// Extracts type
    inline int& Type(int ind){
    	return system->atoms[index[ind]].type;
    }

    /// Extracts residue name
    inline std::string& Resname(int ind){
    	return system->atoms[index[ind]].resname;
    }

    /// Extracts chain
    inline char& Chain(int ind){
    	return system->atoms[index[ind]].chain;
    }

    /// Extracts atom name
    inline std::string& Name(int ind){
    	return system->atoms[index[ind]].name;
    }

    /// Extracts atom mass
    inline float& Mass(int ind){
    	return system->atoms[index[ind]].mass;
    }

    /// Extracts atom charge
    inline float& Charge(int ind){
    	return system->atoms[index[ind]].charge;
    }

    /// Extracts B-factor
    inline float& Beta(int ind){
    	return system->atoms[index[ind]].beta;
    }

    /// Extracts occupancy field
    inline float& Occupancy(int ind){
    	return system->atoms[index[ind]].occupancy;
    }

    /// Extracts residue number
    inline int& Resid(int ind){
    	return system->atoms[index[ind]].resid;
    }

    /// Extracts atom index in the system, which is pointed by selection
    inline int& Index(int ind){
    	return index[ind];
    }

    /// Extracts tag
    inline std::string& Tag(int ind){
    	return system->atoms[index[ind]].tag;
    }

    /// Extracts whole atom
    inline pteros::Atom& Atom(int ind){
        return system->atoms[index[ind]];
    }

    /// Extracts resindex
    inline int& Resindex(int ind){
        return system->atoms[index[ind]].resindex;
    }

    /// Computes VDW radius. Read only.
    inline float VDW(int ind){
        switch(system->atoms[index[ind]].name[0]){
            case 'H': return  0.1;
            case 'C': return  0.17;
            case 'N': return  0.1625;
            case 'O': return  0.149; //mean value used
            case 'S': return  0.1782;
            case 'P': return  0.1871;
            default:  return  0.17;
        }
    }

    /// @}

protected:
    // Indexes of atoms in selection
    std::vector<int> index;
    // Pointer to target system
    System* system;

    // Stores current frame
    int frame;

    // Holds an instance of selection parser
    boost::shared_ptr<Selection_parser> parser;
    void allocate_parser(std::string &sel_text);

    // Private functions for creating selection
    void create_internal(System& sys, std::string& str);
    void create_internal(System& sys, int ind1, int ind2);
    // Private function for deleting selection
    void delete_internal();    

    // Notification responder and connection object
    boost::signals2::connection connection;
    void notify_slot(System_notification code, int b, int e);

    // Here we define read-only accessors for coordinate and mass
    // This is needed because public XYZ() accessor is not const
    // and thus can't be used in const methods
    inline float& _X(int ind) const { return system->traj[frame].coord[index[ind]](0); }
    inline float& _X(int ind, int fr) const { return system->traj[fr].coord[index[ind]](0); }
    inline float& _Y(int ind) const { return system->traj[frame].coord[index[ind]](1); }
    inline float& _Y(int ind, int fr) const { return system->traj[fr].coord[index[ind]](1); }
    inline float& _Z(int ind) const { return system->traj[frame].coord[index[ind]](2); }
    inline float& _Z(int ind, int fr) const { return system->traj[fr].coord[index[ind]](2); }
    inline Eigen::Vector3f& _XYZ(int ind) const { return system->traj[frame].coord[index[ind]]; }
    inline Eigen::Vector3f& _XYZ(int ind, int fr) const { return system->traj[fr].coord[index[ind]]; }
    inline float& _Mass(int ind) const { return system->atoms[index[ind]].mass; }
    inline int _Index(int ind) const { return index[ind]; }
};


}
#endif /* SELECTION_H */
