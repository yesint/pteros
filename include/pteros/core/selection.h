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

#ifndef SELECTION_H
#define SELECTION_H

#include <iostream>

#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "pteros/core/selection_parser.h"
#include "pteros/core/system.h"
#include <boost/shared_ptr.hpp>
#include "pteros/core/typedefs.h"

namespace pteros {

// Forward declaration of friend class
class System;
class Selection_parser;
class Atom_proxy;

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
  // System and Selection are friends because they are closely integrated.
  friend class System;
  friend class Selection_parser;
  friend class Grid_searcher;

  public:
    // Ensure correct 16-bytes-alignment for Eigen sse2 optimizations
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Constructors and basic operators
    /// @{

    /** Main constructor.
    *   @param sys System pointed by this selection
    *   @param str Selection string    
    */
    Selection(System& sys, std::string str);

    /** Constructor with delayed parsing of selection text.
    *   Associates selection with the system @param sys,
        but does not parse selection.
    *   Selection text should be passed later by calling modify().
    */
    Selection(System& sys);

    /// Default constructor for absolutely empty selection.
    Selection();

    /** Constructor, which creates selection from the interval of indexes
        instead of selection string.
        It is much faster then parsing corresponding string, but is limited
        to contigous interval of indexes.
        @param sys System pointed by this selection
        @param ind1 First index in interval
        @param ind2 Last index in interval (inclusive!)
     */
    Selection(System& sys, int ind1, int ind2);

    /// Assignment operator
    Selection& operator=(Selection sel);

    /// Copy constructor
    Selection(const Selection& sel);

    /// Equality operator
    /// Selection are compared by their indexes
    bool operator==(const Selection &other) const;

    /// Inequality operator
    /// /// Selection are compared by their indexes
    bool operator!=(const Selection &other) const {
        return !(*this == other);
    }

    /// Indexing operator. Returns an Atom_proxy object,
    /// which incapsulates atom and its coordinates for current frame
    /// Can only be used as r-value (can not be assigned to).
    Atom_proxy operator[](int ind) const;

    /// Destructor
    ~Selection();

    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Modification of existing selection
    /// @{

    /// Append another selection to this one. Acts like logical "OR" between selections.
    void append(Selection& sel);

    /// Append absolute index to selection
    void append(int ind);

    /// Modifies both system and string in selection.    
    void modify(System& sys, std::string str);

    /** Modifies selection string in existing selection.
    *   @param str New value of selection text. Selection is re-parsed immediately with
    *   this new value.
    */
    void modify(std::string str);

    /// Change system for selection (clears selection)
    void modify(System& sys);

    /// Modifies both system and selection using the range of indexes.   
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

    /** Clears selection and frees memory, but do not delete it.
    *   Selection remains registered in parent system, but is cleared from any data.
    *   Subsequent call of modify() may populate it again.
    */
    void clear();
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Iterator access
    /// @{

    /// Get const iterator for begin of index
    std::vector<int>::const_iterator index_begin() const { return index.begin(); }

    /// Get const iterator for the end of index
    std::vector<int>::const_iterator index_end() const { return index.end(); }

    /// Get back_insert_iterator for index
    std::back_insert_iterator<std::vector<int> > index_back_inserter() {
        return std::back_inserter(index);
    }

    /// Forward random-access iterator for Selection class.
    /// Iterator reffers to Atom_proxy object
    class iterator;

    /// Begin iterator
    iterator begin();
    /// End iterator
    iterator end();
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Get and set functions
    /// @{

    /// Get current frame selection is pointing to
    int get_frame() const {return frame;}

    /// Set current frame for selection
    void set_frame(int fr);

    /// Get pointer to the system, which owns this selection
    System* get_system() const { return system; }

    /// Get selection text
    std::string get_text() const;

    /// Get vector of all indexes in selection
    std::vector<int> get_index() const { return index; }    

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

    /// Get coordinates of all atoms in this selection for the current frame
    Eigen::MatrixXf get_xyz() const;

    /// Get coordinates of all atoms in this selection for the current frame
    void get_xyz(MatrixXf_ref res) const;

    /// Set coordinates of this selection for current frame
    void set_xyz(MatrixXf_const_ref coord);

    /// Get masses of all atoms in selection
    std::vector<float> get_mass() const;

    /// Set atom masses in selection to the values from supplied vector.
    /// Its size must be the save as the size of selection.
    void set_mass(const std::vector<float> m);

    /// Sets masses of all selected atoms to the same given value.
    void set_mass(float data);        

    /// Get beta
    std::vector<float> get_beta() const;

    /// Set beta in selection to the values from supplied vector.
    /// Its size must be the save as the size of selection.
    void set_beta(std::vector<float>& data);

    /// Sets beta of all selected atoms to the same given value.
    void set_beta(float data);    
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Computing properties of selection
    /// @{

    /** Get the center of selection.
    * @param mass_weighted Use mass-weighting
    * @param periodic Account for periodic boundary conditions.
    * Please note that if the size of selection is larger than 1/2 of the box size in
    * any dimension you will get incorrect results if periodic is set to true.
    */
    Eigen::Vector3f center(bool mass_weighted = false, bool periodic = false) const;

    /// Get minimal and maximal coordinates in selection
    void minmax(Vector3f_ref min, Vector3f_ref max) const;

#ifdef USE_POWERSASA
    /// Get the SASA. Easy way - only returns SASA area of selection
    float sasa(float probe_r = 0.14) const;

    /// Get the SASA. Detailed way - returns area and computes volume and per-atom values
    float sasa(float probe_r = 0.14, float* total_volume = NULL,
               std::vector<float>* area_per_atom = NULL,
               std::vector<float>* volume_per_atom = NULL) const;
#endif

    /// Computes average structure over the range of frames
    Eigen::MatrixXf average_structure(int b=0, int e=-1) const;

    /** Extracts X,Y,Z for given atom index for specified range of frames
        (gets trajectory of given atom).
    *   Result is returned as MatrixXf, where i-th column is an XYZ vector for frame i.
    */
    Eigen::MatrixXf atom_traj(int ind, int b=0, int e=-1) const;

    /// Computes the central momens of inertia and principal axes of inertia
    void inertia(Vector3f_ref moments, Matrix3f_ref axes, bool periodic = false) const;

    /// Computes radius of gyration for selection
    float gyration(bool periodic = false) const;
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Geometry transformation functions
    /// @{

    /// Translate selection by given vector
    void translate(Vector3f_const_ref v);

    /// Rotate selection around given axis relative to center of masses
    /// @param axis Axis of rotation 0=X, 1=Y, 2=Z
    /// @param angle Rotation angle in radians
    void rotate(int axis, float angle);

    /// Rotate selection around given axis relative to given pivot
    /// @param axis Axis of rotation 0=X, 1=Y, 2=Z
    /// @param angle Rotation angle in radians
    /// @param pivot Rotation around this pivot
    void rotate(int axis, float angle, Vector3f_const_ref pivot);

    /// Rotate selection around given vector relative to given pivot
    /// @param direction Rotate around this vector
    /// @param angle Rotation angle in radians
    /// @param pivot Rotation around this pivot
    void rotate(Vector3f_const_ref direction, float angle, Vector3f_const_ref pivot);

    /// Rotation with the given 3x3 rotation matrix around point (0,0,0)
    void rotate(Matrix3f_const_ref m);

    /// Rotation by given angles around X, Y and Z with given pivot
    void rotate(Vector3f_const_ref angles, Vector3f_const_ref pivot);

    /// Wraps whole selection to the periodic box
    void wrap(Vector3i_const_ref dims = Eigen::Vector3i::Ones());

    /** Unwraps selection to make it whole if possible (without jumps over periodic box boundary).
     * The periodic center of mass is used as an anchor point.
     * Please note that if the size of selection is larger than 1/2 of the box size in
     * any dimension unwrap() will not work as expected and will not make selection "compact"!
    */
    void unwrap(Vector3i_const_ref dims = Eigen::Vector3i::Ones());

    /** Unwraps selection to make it whole (without jumps over periodic box boundary).
     * based on preserving all bonds. The maximal bond length is given by d.
     * This method works reliably in any case, but is much slower than unwrap()
     */
    void unwrap_bonds(float d = 0.2, Vector3i_const_ref dims = Eigen::Vector3i::Ones());

    /** Get transform for orienting selection by principal axes.
     * Please note that if the size of selection is larger than 1/2 of the box size in
     * any dimension you will get funny results if is_periodic is set to true.
     */    
    Eigen::Affine3f principal_transform(bool is_periodic = false) const;

    /** Orient molecule by its principal axes.
     * The same as
     * \code
     * Eigen::Affine3f tr = sel.principal_transform();
     * sel.apply_transform(tr);
     * \endcode
     */
    void principal_orient(bool is_periodic = false);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Fitting and RMSD functions
    /// @{

    /// RMSD between current and another frame
    float rmsd(int fr) const;

    /// RMSD between two frames
    float rmsd(int fr1, int fr2) const;

    /// RMSD between two selections of the same size
    friend float rmsd(const Selection& sel1, const Selection& sel2);

    /** RMSD between two selections of the same size (for given frames)
    *   @param sel1 First selection
    *   @param fr1 Frame for first selection
    *   @param sel2 Second selection
    *   @param fr2 Frame for second selection
    */
    friend float rmsd(const Selection& sel1, int fr1, const Selection& sel2, int fr2);

    /// Fit two selection of the same size. sel1 is modified to be fit to sel2.
    friend void fit(Selection& sel1, const Selection& sel2);

    /// Fit all frames in the trajectory to reference frame
    void fit_trajectory(int ref_frame=0, int b=0, int e=-1);

    /// Returns fitting transformation for two given selections of the same size
    friend Eigen::Affine3f fit_transform(const Selection& sel1, const Selection& sel2);

    /// Returns fit transformation between frames fr1 and fr2
    Eigen::Affine3f fit_transform(int fr1, int fr2) const;

    /// Fits frame fr1 to fr2
    void fit(int fr1, int fr2);

    /// Apply fitting transformation
    void apply_transform(Eigen::Affine3f& t);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Energy functions
    /// @{

    /// Self-energy of selection computed within given interaction cut-off.
    /// If cutoff is 0 or negative computes interaction of all atoms
    /// (very slow for large selections)
    Energy_components non_bond_energy(float cutoff=0.25, bool periodic=true) const;

    /// Non-bond energy between two selections computed within given interaction cut-off.
    /// If cutoff is 0 or negative computes all pairs of atoms (very slow for large selections)
    /// fr = -1 computes for current frame of selection 1.
    friend Energy_components non_bond_energy(const Selection& sel1,
                                             const Selection& sel2,
                                             float cutoff = 0.25,
                                             int fr = -1,
                                             bool periodic = true);

    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name File IO
    /// @{

    /** Write structure or trajectory for selection. Type is determined by extension.
    *   Frames from b to e are written.
    *   If b and e are not set they default to current frame
    */
    // Can't be made const because of internal calls
    void write(std::string fname,int b=-1,int e=-1);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Building functions
    /// @{

    /// Duplicate current selection in the parent system
    /// @param res_sel Pointer to selection of duplicated atoms
    void atoms_dup(Selection* res_sel = NULL);

    /// Delete all atoms of current selection from the parent system
    /// Current selection becomes invalid after this!
    void atoms_delete();

    /// Creates multiple copies of selection in the parent system and
    /// distributes them in a grid
    void distribute(Vector3i_const_ref ncopies, Vector3f_const_ref shift);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Utility functions
    /// @{

    /// Get the size of selection
    int size() const {return index.size();}

    /// Returnss true if selection was created from text string and false if it was
    /// constructed 'by hand' by appending indexes or other selections
    bool text_based() const {
        return sel_text!="";
    }

    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Splitting selection into parts
    /// @{

    /// Split current selection into several selections according to
    /// the interatomic distances. Each resulting selection is a group
    /// of atoms connected by distances less than d
    void split_by_connectivity(float d, std::vector<Selection>& res);

    /// Split selection by residue index
    void split_by_residue(std::vector<Selection>& res);

    /// Selects each residue, which is referenced by selection.
    /// All selections for residues are placed into supplied vector.
    /// Handles multiple chains correctly.    
    void each_residue(std::vector<Selection>& sel) const;
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** @name Inline accessors
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

    inline const float& X(int ind) const {
        return system->traj[frame].coord[index[ind]](0);
    }

    /// Extracts X for given frame frame fr
    inline float& X(int ind, int fr){
        return system->traj[fr].coord[index[ind]](0);
    }

    inline const float& X(int ind, int fr) const {
        return system->traj[fr].coord[index[ind]](0);
    }

    /// Extracts Y for current frame
    inline float& Y(int ind){
        return system->traj[frame].coord[index[ind]](1);
    }

    inline const float& Y(int ind) const {
        return system->traj[frame].coord[index[ind]](1);
    }

    /// Extracts Y for given frame frame fr
    inline float& Y(int ind, int fr){
    	return system->traj[fr].coord[index[ind]](1);
    }

    inline const float& Y(int ind, int fr) const {
        return system->traj[fr].coord[index[ind]](1);
    }

    /// Extracts Z for current frame
    inline float& Z(int ind){
    	return system->traj[frame].coord[index[ind]](2);
    }

    inline const float& Z(int ind) const {
        return system->traj[frame].coord[index[ind]](2);
    }

    /// Extracts Z for given frame frame fr
    inline float& Z(int ind, int fr){
    	return system->traj[fr].coord[index[ind]](2);
    }

    inline const float& Z(int ind, int fr) const {
        return system->traj[fr].coord[index[ind]](2);
    }

    /// Extracts X,Y and Z for current frame
    inline Eigen::Vector3f& XYZ(int ind){
    	return system->traj[frame].coord[index[ind]];
    }

    inline const Eigen::Vector3f& XYZ(int ind) const {
        return system->traj[frame].coord[index[ind]];
    }

    /// Extracts X,Y and Z for given frame frame fr
    inline Eigen::Vector3f& XYZ(int ind, int fr){
    	return system->traj[fr].coord[index[ind]];
    }

    inline const Eigen::Vector3f& XYZ(int ind, int fr) const {
        return system->traj[fr].coord[index[ind]];
    }

    /// Extracts type
    inline int& Type(int ind){
    	return system->atoms[index[ind]].type;
    }

    inline const int& Type(int ind) const {
        return system->atoms[index[ind]].type;
    }

    /// Extracts typename
    inline std::string& Type_name(int ind){
        return system->atoms[index[ind]].type_name;
    }

    inline const std::string& Type_name(int ind) const {
        return system->atoms[index[ind]].type_name;
    }

    /// Extracts residue name
    inline std::string& Resname(int ind){
    	return system->atoms[index[ind]].resname;
    }

    inline const std::string& Resname(int ind) const {
        return system->atoms[index[ind]].resname;
    }

    /// Extracts chain
    inline char& Chain(int ind){
    	return system->atoms[index[ind]].chain;
    }

    inline const char& Chain(int ind) const {
        return system->atoms[index[ind]].chain;
    }

    /// Extracts atom name
    inline std::string& Name(int ind){
    	return system->atoms[index[ind]].name;
    }

    inline const std::string& Name(int ind) const {
        return system->atoms[index[ind]].name;
    }

    /// Extracts atom mass
    inline float& Mass(int ind){
    	return system->atoms[index[ind]].mass;
    }

    inline const float& Mass(int ind) const {
        return system->atoms[index[ind]].mass;
    }

    /// Extracts atom charge
    inline float& Charge(int ind){
    	return system->atoms[index[ind]].charge;
    }

    inline const float& Charge(int ind) const {
        return system->atoms[index[ind]].charge;
    }

    /// Extracts B-factor
    inline float& Beta(int ind){
    	return system->atoms[index[ind]].beta;
    }

    inline const float& Beta(int ind) const {
        return system->atoms[index[ind]].beta;
    }

    /// Extracts occupancy field
    inline float& Occupancy(int ind){
    	return system->atoms[index[ind]].occupancy;
    }

    inline const float& Occupancy(int ind) const {
        return system->atoms[index[ind]].occupancy;
    }

    /// Extracts residue number
    inline int& Resid(int ind){
    	return system->atoms[index[ind]].resid;
    }

    inline const int& Resid(int ind) const {
        return system->atoms[index[ind]].resid;
    }

    /// Extracts atom index in the system, which is pointed by selection
    inline int& Index(int ind){
    	return index[ind];
    }

    inline const int& Index(int ind) const {
        return index[ind];
    }

    /// Extracts tag
    inline std::string& Tag(int ind){
    	return system->atoms[index[ind]].tag;
    }

    inline const std::string& Tag(int ind) const {
        return system->atoms[index[ind]].tag;
    }

    /// Extracts whole atom
    inline pteros::Atom& Atom_data(int ind){
        return system->atoms[index[ind]];
    }

    inline const pteros::Atom& Atom_data(int ind) const {
        return system->atoms[index[ind]];
    }

    /// Extracts resindex
    inline int& Resindex(int ind){
        return system->atoms[index[ind]].resindex;
    }

    inline const int& Resindex(int ind) const {
        return system->atoms[index[ind]].resindex;
    }

    /// Computes VDW radius. Read only.
    inline float VDW(int ind) const {
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
    // Row text of selection
    std::string sel_text;    
    // Indexes of atoms in selection
    std::vector<int> index;
    // Pointer to target system
    System* system;

    // Stores current frame
    int frame;

    // Holds an instance of selection parser
    boost::shared_ptr<Selection_parser> parser;
    void allocate_parser();

    // Private functions for creating selection
    void create_internal(System& sys, std::string& str);
    void create_internal(System& sys, int ind1, int ind2);
    // Private function for deleting selection
    void delete_internal();        
};

//==============================================================================

/// Auxilary type used to incapsulate the atom and its current coordinates
/// Used internally in Selection::operator[] and in iterator access to Selection
/// Objects of this class should not be created by the user in normal situation.
class Atom_proxy {
    friend class Selection::iterator;
public:
    Atom_proxy(){}
    Atom_proxy(Selection* s, int i): sel(s), ind(i) {}        

    /// Accessors. Const and non-const versions.
    inline int& Resid(){ return sel->Resid(ind); }
    inline const int& Resid() const { return sel->Resid(ind); }

    inline int& Index(){ return sel->Index(ind); }
    inline const int& Index() const { return sel->Index(ind); }

    inline std::string& Name(){ return sel->Name(ind); }
    inline const std::string& Name() const { return sel->Name(ind); }

    inline char& Chain(){ return sel->Chain(ind); }
    inline const char& Chain() const { return sel->Chain(ind); }

    inline std::string& Resname(){ return sel->Resname(ind); }
    inline const std::string& Resname() const { return sel->Resname(ind); }

    inline std::string& Tag(){ return sel->Tag(ind); }
    inline const std::string& Tag() const { return sel->Tag(ind); }

    inline float& Occupancy(){ return sel->Occupancy(ind); }
    inline const float& Occupancy() const { return sel->Occupancy(ind); }

    inline float& Beta(){ return sel->Beta(ind); }
    inline const float& Beta() const { return sel->Beta(ind); }

    inline int& Resindex(){ return sel->Resindex(ind); }
    inline const int& Resindex() const { return sel->Resindex(ind); }

    inline float& Mass(){ return sel->Mass(ind); }
    inline const float& Mass() const { return sel->Mass(ind); }

    inline float& Charge(){ return sel->Charge(ind); }
    inline const float& Charge() const { return sel->Charge(ind); }

    inline int& Type(){ return sel->Type(ind); }
    inline const int& Type() const { return sel->Type(ind); }

    inline std::string& Type_name(){ return sel->Type_name(ind); }
    inline const std::string& Type_name() const { return sel->Type_name(ind); }

    inline float& X(){ return sel->X(ind); }
    inline const float& X() const { return sel->X(ind); }

    inline float& Y(){ return sel->Y(ind); }
    inline const float& Y() const { return sel->Y(ind); }

    inline float& Z(){ return sel->Z(ind); }
    inline const float& Z() const { return sel->Z(ind); }

    inline Eigen::Vector3f& XYZ(){ return sel->XYZ(ind); }
    inline const Eigen::Vector3f& XYZ() const { return sel->XYZ(ind); }

    inline float& X(int fr){ return sel->X(ind,fr); }
    inline const float& X(int fr) const { return sel->X(ind,fr); }

    inline float& Y(int fr){ return sel->Y(ind,fr); }
    inline const float& Y(int fr) const { return sel->Y(ind,fr); }

    inline float& Z(int fr){ return sel->Z(ind,fr); }
    inline const float& Z(int fr) const { return sel->Z(ind,fr); }

    inline Eigen::Vector3f& XYZ(int fr){ return sel->XYZ(ind,fr); }
    inline const Eigen::Vector3f& XYZ(int fr) const { return sel->XYZ(ind,fr); }

    inline Atom& Atom_data(){ return sel->Atom_data(ind); }
    inline const Atom& Atom_data() const { return sel->Atom_data(ind); }

    /// Equality operator
    bool operator==(const Atom_proxy& other) const {
        return (sel==other.sel && ind==other.ind);
    }

    /// Inequality operator
    bool operator!=(const Atom_proxy &other) const {
        return !(*this == other);
    }

protected:
    Selection* sel;
    int ind;
};

//==============================================================================

/// Random-access forward iterator for Selection
class Selection::iterator {
public:
    typedef Atom_proxy value_type;
    typedef int difference_type;
    typedef Atom_proxy* pointer;
    typedef Atom_proxy& reference;
    typedef std::forward_iterator_tag iterator_category;

    iterator(Selection* sel, int pos) { proxy.sel = sel; proxy.ind = pos; }
    iterator operator++() { iterator tmp = *this; proxy.ind++; return tmp; }
    iterator operator++(int junk) { proxy.ind++; return *this; }
    Atom_proxy& operator*() { return proxy; }
    Atom_proxy* operator->() { return &proxy; }
    bool operator==(const iterator& rhs) { return proxy == rhs.proxy; }
    bool operator!=(const iterator& rhs) { return proxy != rhs.proxy; }
private:
    Atom_proxy proxy;
};

} // namespace pteros
#endif /* SELECTION_H */
