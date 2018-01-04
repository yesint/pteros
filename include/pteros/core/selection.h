/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2015, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite:
 *
 *     Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source
 *     C++ library for molecular analysis",
 *     Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/

#ifndef SELECTION_H
#define SELECTION_H

#include <iostream>

#include <string>
#include <memory>
#include <map>
#include <functional>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "pteros/core/system.h"
#include "pteros/core/typedefs.h"

namespace pteros {

// Forward declaration of friend classes
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
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Constructors and operators
    /// @{

    /// Default constructor for absolutely empty selection.
    Selection();

    /** Associates selection with the system @param sys,
    *   Selection content should be set later by calling modify().
    */
    Selection(const System& sys);

    /** Main constructor.
    *   @param sys System pointed by this selection
    *   @param str Selection string    
    */
    Selection(const System& sys, std::string str, int fr = 0);

    /** Constructor, which creates selection from the interval of indexes
        instead of selection string.
        It is much faster then parsing corresponding string, but is limited
        to contigous interval of indexes.
        @param sys System pointed by this selection
        @param ind1 First index in interval
        @param ind2 Last index in interval (inclusive!)        
     */
    Selection(const System& sys, int ind1, int ind2);

    /// Constructor from vector of indexes
    /// Vector may be in any order and may contain duplicates.
    Selection(const System& sys, const std::vector<int>& ind);

    /// Constructor from the pair of iterators to int sequence
    /// Sequence may be in any order and may contain duplicates.
    Selection(const System& sys,
              std::vector<int>::iterator it1,
              std::vector<int>::iterator it2);

    /** Constructor which takes user-defined callback
      Callback takes the system as first argument, target frame number as the second
      and the vector to be filled by selected atom indexes.
      Vector may be filled in any order and may contain duplicates.
      \warning
      Resulting selection is neither coordinate-dependent nor text-based.
      It will not recompute itself on the frame change even if it involves atom coordinates.
    */
    Selection(const System& sys,
              const std::function<void(const System&,int,std::vector<int>&)>& callback,
              int fr = 0);

    /// Copy constructor
    Selection(const Selection& sel);

    /// Destructor
    virtual ~Selection();

    /// Assignment operator
    Selection& operator=(Selection sel);

    /// Equality operator
    /// Selection are compared by their indexes
    bool operator==(const Selection &other) const;

    /// Inequality operator
    /// Selection are compared by their indexes
    bool operator!=(const Selection &other) const {
        return !(*this == other);
    }

    /** Indexing operator. Returns an Atom_proxy object,
     which incapsulates atom and its coordinates for current frame
     Can only be used as r-value (can not be assigned to)
     however individual fields of the proxy object are writable

     \code
     Selection sel(sys,"name CA");
     sel[1] = sel[0]; // ERROR! Can't assign to proxy object returned by []!
     sel[0].Name() = "A"; // OK to write fields of proxy object
     \endcode

     Main goal of [] operator is usage in range-based for loops:
     \code
     Selection sel(sys,"name CA");
     for(auto a: sel){
        cout << a.Name() << " " << a.X() << endl;
     }
     \endcode

     Otherwise it is slower than conventional syntax like sel.Name(i)
     */
    Atom_proxy operator[](int ind) const;

    /// Writing selection to stream.
    /// Outputs indexes as a space separated list
    friend std::ostream& operator<<(std::ostream& os, const Selection& sel);

    /// Creates new Selection, which is the logical OR of two parent selections.
    /// Parent selections are not modified.    
    friend Selection operator|(const Selection& sel1, const Selection& sel2);        

    /// Creates new Selection, which is the logical AND of two parent selections.
    /// Parent selections are not modified.
    friend Selection operator&(const Selection& sel1, const Selection& sel2);

    /// Creates new Selection, by removing all atoms of sel2 from sel1.
    /// Parent selections are not modified.
    /// \warning This operator is \em not commutative!
    friend Selection operator-(const Selection& sel1, const Selection& sel2);

    /// Creates new Selection, which is a logical negation of existing one.
    /// Parent selection is not modified
    Selection operator~() const;        
    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Modification of existing selection
    /// @{

    /// Append another selection to this one.
    void append(const Selection& sel);

    /// Append absolute index to selection
    void append(int ind);

    /// Remove all atoms of sel from current selection
    void remove(const Selection& sel);

    /// Remove given absolute index from current selection
    void remove(int ind);

    /// Inverts selection in place by selection those atoms which were not selected
    void invert();

    /// Sets new system for selection.
    /// \warning This clears selection index and leaves it empty!
    void set_system(const System& sys);

    /** Modifies selection string in existing selection.
    *   @param str New value of selection text. Selection is re-parsed immediately with
    *   this new value.
    */
    void modify(std::string str, int fr = 0);

    /// Modifies selection using the range of indexes
    void modify(int ind1, int ind2);

    /// Modifies selection using vector of indexes
    /// Vector may be in any order and may contain duplicates.
    void modify(const std::vector<int>& ind);

    /// Modifies selection using pair of iterators to index vector
    /// Vector may be in any order and may contain duplicates.
    void modify(std::vector<int>::iterator it1, std::vector<int>::iterator it2);

    /** Modifies selection using user-defined callback.
      Callback takes the system as first argument, target frame number as the second
      and the vector to be filled by selected atom indexes.
      Vector may be filled in any order and may contain duplicates.
      \warning
      Resulting selection is neither coordinate-dependent nor text-based.
      It will not recompute itself on the frame change even if it involves atom coordinates.
    */
    void modify(const std::function<void(const System&,int,std::vector<int>&)>& callback, int fr = 0);

    /// Convenience function, which combines set_system and modify(str)
    void modify(const System& sys, std::string str, int fr = 0);

    /// Convenience function, which combines set_system and modify(int,int)
    void modify(const System& sys, int ind1, int ind2);

    /// Convenience function, which combines set_system and modify(vector<int>)
    void modify(const System& sys, const std::vector<int>& ind);

    /// Convenience function, which combines set_system and modify(iter,iter)
    void modify(const System& sys,
                std::vector<int>::iterator it1,
                std::vector<int>::iterator it2);

    /// Convenience function, which combines set_system and modify(callback)
    void modify(const System& sys, const std::function<void(const System&,int,std::vector<int>&)>& callback, int fr = 0);

    /** Recomputes selection without re-parsing selection text.
    *   Only makes sense for coordinate-dependent selections when the coordinates change.
    *   Called authomatically for coordinate-dependent selections by set_frame()
    *   If selection is not coordinate-dependent does nothing.
    */
    void apply();

    /** Recomputes selection completely.
    *   May be used when new file is loaded into the system, or when atoms are
    *   created/deleted. Forces re-parsing of selection text.
    *   For non-textual selections does nothing.
    */
    void update();

    /** Clears selection and frees memory.
    *   Selection is still assigned to its system after clearing.
    *   Subsequent call of modify() may populate it again.
    */
    void clear();
    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Sub-selections
    /**
      Sub-selections allow selecting atoms \e inside existing selection
      (narrowing or refining existing selection in other terms).
      Sub-selections could be very useful in the following situation.
      Suppose that we need to create separate selections for N,C,CA and CB atoms
      of particular protein residue. With "normal" selections the following code could be used:
      \code
      Selection sel_N(sys,"protein and resid 1 and name N");
      Selection sel_C(sys,"protein and resid 1 and name C");
      Selection sel_CA(sys,"protein and resid 1 and name CA");
      Selection sel_CB(sys,"protein and resid 1 and name CB");
      \endcode
      The problem with this code is that we are looping over \e all atoms in the system four times,
      ones in each selection. This is very inefficient since we only need to find our
      residue with "protein and resid 1" (one loop over all atoms) and then we need to search
      \e inside this residue four times (looping over ~10 atoms only). This problem is not
      apparent for small systems but becomes very painful for the systems with millions of atoms.
      Subselections solve this problem:
      \code
      Selection residue1(sys,"protein and resid 1");
      auto sel_N = residue1.select("name N");
      auto sel_C = residue1.select("name C");
      auto sel_CA = residue1.select("name CA");
      auto sel_CB = residue1.select("name CB");
      \endcode

      Subselections inherit the system and frame from the parent. The search in sub-selections
      is performed over selected atoms of the parent only (the only exception from this rule are
      within selections which involve seacrh over all atoms by design).
    */
    /// @{
    Selection select(std::string str);
    Selection operator()(std::string str);

    /// Local selection indexes are used!
    Selection select(int ind1, int ind2);
    Selection operator()(int ind1, int ind2);

    /// Local selection indexes are used!
    Selection select(const std::vector<int>& ind);
    Selection operator()(const std::vector<int>& ind);
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
    /// \note If selection is coordinate-dependent it is re-evaluated by calling apply()
    void set_frame(int fr);

    /// Get pointer to the system, which is pointed by this selection
    System* get_system() const { return system; }

    /// Get selection text
    /// If selection is not textual returns generated index-based text "index i j k..."
    std::string get_text() const;

    /// Get vector of all indexes in selection
    std::vector<int> get_index() const { return index; }    

    /// Get vector of all chains in selection
    std::vector<char> get_chain(bool unique=false) const;

    /// Set chains from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_chain(const std::vector<char>& data);

    /// Sets chain of all selected atoms to the same given value.
    void set_chain(char data);


    /// Get vector of all resid's in selection
    /// \warning Same resid's could be present in different chains!
    std::vector<int> get_resid(bool unique=false) const;

    /// Set resid's in selection from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_resid(const std::vector<int>& data);

    /// Sets resid of all selected atoms to the same given value.
    void set_resid(int data);


    /// Get vector of all atom names in selection
    std::vector<std::string> get_name(bool unique=false) const;

    /// Set atom names in selection from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_name(const std::vector<std::string>& data);

    /// Sets atom names of all selected atoms to the same given value.
    void set_name(std::string data);


    /// Get vector of all resnames in selection
    std::vector<std::string> get_resname(bool unique=false) const;

    /// Set resnames in selection from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_resname(const std::vector<std::string>& data);

    /// Sets resnames of all selected atoms to the same given value.
    void set_resname(std::string data);


    /// Get vector of all resindexes in selection. Resindexes are unique
    /// regardless the number of the chains.
    std::vector<int> get_resindex(bool unique=false) const;


    /// Get coordinates of all atoms in this selection for the current frame
    /// By default columnt of the matrix contain atom coordinates.
    /// Optional flag turns this to row-major format (mainly used by python binding)
    Eigen::MatrixXf get_xyz(bool make_row_major_matrix = false) const;

    /// Get coordinates of all atoms in this selection for the current frame in existing matrix
    void get_xyz(MatrixXf_ref res) const;

    /// Set coordinates of this selection for current frame
    void set_xyz(MatrixXf_const_ref coord);


    /// Get masses of all atoms in selection
    std::vector<float> get_mass() const;

    /// Set atom masses in selection to the values from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_mass(const std::vector<float>& m);

    /// Sets masses of all selected atoms to the same given value.
    void set_mass(float data);        


    /// Get beta
    std::vector<float> get_beta() const;

    /// Set beta in selection to the values from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_beta(const std::vector<float>& data);

    /// Sets beta of all selected atoms to the same given value.
    void set_beta(float data);


    /// Get occupancy
    std::vector<float> get_occupancy() const;

    /// Set occupancy in selection to the values from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_occupancy(const std::vector<float>& data);

    /// Sets occupancy of all selected atoms to the same given value.
    void set_occupancy(float data);


    /// Get tags
    std::vector<std::string> get_tag(bool unique=false) const;

    /// Set tags in selection to the values from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_tag(const std::vector<std::string>& data);

    /// Set tags of all selected atoms to the same given value.
    void set_tag(std::string data);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Computing properties of selection
    /// @{

    /** Get the center of selection.
     @param mass_weighted Use mass-weighting
     @param periodic Account for periodic boundary conditions.

     \warning
     If the size of selection is larger than 1/2 of the box size in
     any dimension you will get incorrect results if periodic is set to true.
     This is not checked automatically!
     In this case use one of unwrapping options first.
    */
    Eigen::Vector3f center(bool mass_weighted = false, bool periodic = false, int leading_index = 0) const;

    /// Get minimal and maximal coordinates in selection
    void minmax(Vector3f_ref min, Vector3f_ref max) const;

    /// Get the SASA using powersasa algorithm. Returns area and computes volume and per-atom values if asked
    float powersasa(float probe_r = 0.14,
               std::vector<float>* area_per_atom = nullptr,
               float* total_volume = nullptr,
               std::vector<float>* volume_per_atom = nullptr) const;

    /// Get SASA using Shrake and Rupley algorithm (slower than powersasa)
    float sasa(float probe_r = 0.14, std::vector<float>* area_per_atom = nullptr, int n_sphere_points = 960) const;

    /// Computes average structure over the range of frames
    Eigen::MatrixXf average_structure(int b=0, int e=-1, bool make_row_major_matrix = false) const;

    /** Extracts X,Y,Z for given atom index for specified range of frames
        (gets trajectory of given atom).
    *   Result is returned as MatrixXf, where i-th column (or row) is an XYZ vector for frame i.
    */
    Eigen::MatrixXf atom_traj(int ind, int b=0, int e=-1, bool make_row_major_matrix = false) const;

    /** Computes the central momens of inertia and principal axes of inertia
     \warning
     If the size of selection is larger than 1/2 of the box size in
     any dimension you will get incorrect results if periodic is set to true.
     This is not checked automatically!
     In this case use one of unwrapping options first.
    */
    void inertia(Vector3f_ref moments, Matrix3f_ref axes, bool periodic = false, bool leading_index = 0) const;

    /** Computes radius of gyration for selection
     \warning
     If the size of selection is larger than 1/2 of the box size in
     any dimension you will get incorrect results if periodic is set to true.
     This is not checked automatically!
     In this case use one of unwrapping options first.
    */
    float gyration(bool periodic = false, bool leading_index = 0) const;

    /// Get distance between two atoms (periodic in given dimensions if needed).
    /// \note
    /// This function takes selection indexes, not absolute indexes.
    float distance(int i, int j, bool is_periodic = true,
                   Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// Get angle in degrees between three atoms (periodic in given dimensions if needed).
    /// \note
    /// This function takes selection indexes, not absolute indexes.
    float angle(int i, int j, int k, bool is_periodic = true,
                Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

    /// Get dihedral angle in degrees between three atoms (periodic in given dimensions if needed).
    /// \note
    /// This function takes selection indexes, not absolute indexes.
    float dihedral(int i, int j, int k, int l, bool is_periodic = true,
                Vector3i_const_ref dims = Eigen::Vector3i::Ones()) const;

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
      The periodic center of mass is used as an anchor point.
      \warning
      If the size of selection is larger than 1/2 of the box size in
      any dimension unwrap() will not work as expected and will not make selection "compact"!
    */
    void unwrap(int leading_index = -1, Vector3i_const_ref dims = Eigen::Vector3i::Ones());

    /** Unwraps selection to make it whole (without jumps over periodic box boundary).
     * based on preserving all bonds.
     * This method works reliably in any case, but is much slower than unwrap()
     * @param d Maximal bond length.
     * @param leading_index Local index of the reference atom, which doesn't move.
     * @return Number of disconnected pieces after unwrapping. 1 means solid selection.
     */
    int unwrap_bonds(float d = 0.2, int leading_index = 0,
                     Vector3i_const_ref dims = Eigen::Vector3i::Ones());

    /** Get transform for orienting selection by principal axes.
     * \warning
     * If the size of selection is larger than 1/2 of the box size in
     * any dimension you will get funny results if @param is_periodic is set to true.
     */    
    Eigen::Affine3f principal_transform(bool is_periodic = false, bool leading_index = 0) const;

    /** Orient molecule by its principal axes.
     * The same as
     * \code
     * Eigen::Affine3f tr = sel.principal_transform();
     * sel.apply_transform(tr);
     * \endcode
     */
    void principal_orient(bool is_periodic = false, bool leading_index = 0);
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
    *   @param fr1 %Frame for first selection
    *   @param sel2 Second selection
    *   @param fr2 %Frame for second selection
    */
    friend float rmsd(const Selection& sel1, int fr1, const Selection& sel2, int fr2);

    /// Fit two selection of the same size. sel1 is modified to be fit to sel2.
    friend void fit(Selection& sel1, const Selection& sel2);

    /// Fit specified frames in the trajectory to reference frame
    void fit_trajectory(int ref_frame=0, int b=0, int e=-1);

    /// Returns fitting transformation for two given selections of the same size
    friend Eigen::Affine3f fit_transform(const Selection& sel1, const Selection& sel2);

    /// Returns fit transformation between frames fr1 and fr2
    Eigen::Affine3f fit_transform(int fr1, int fr2) const;

    /// Fits frame fr1 to fr2
    void fit(int fr1, int fr2);

    /// Apply fitting transformation
    void apply_transform(const Eigen::Affine3f& t);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Energy functions
    /// @{

    /// Self-energy of selection computed within given interaction cut-off.
    /// If cutoff is 0 the cutoff from topology is used.
    Eigen::Vector2f non_bond_energy(float cutoff=0, bool periodic=true) const;

    /// Non-bond energy between two selections computed within given interaction cut-off.
    /// If cutoff is 0 the cutoff from topology is used.
    /// fr = -1 computes for current frame of selection 1.
    friend Eigen::Vector2f non_bond_energy(const Selection& sel1,
                                             const Selection& sel2,
                                             float cutoff = 0,
                                             int fr = -1,
                                             bool periodic = true);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name File IO
    /// @{

    /** Write structure or trajectory for selection. File type is determined by extension.
    *   Frames from b to e are written.
    *   If @param b is not set or -1 it means current frame
    *   If @param e is not set or -1 it means the last frame
    */
    // Can't be made const because of internal calls
    void write(std::string fname, int b=-1,int e=-1);

    void write(const std::unique_ptr<Mol_file>& handler, Mol_file_content what,int b=-1,int e=-1);    
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Utility functions
    /// @{

    /// Get the size of selection
    int size() const {return index.size();}

    /// Returns true if selection was created from text string and false if it was
    /// constructed 'by hand' by appending indexes or other selections
    bool text_based() const {
        return sel_text!="";
    }

    /// Returns true if selection is coordinate-dependent and is able to recompute
    /// itself on the change of frame
    bool coord_dependent() const {
        return (bool)parser;
    }

    /// "Flattens" selection by removing coordinate dependence and making it not text-based.
    /// Resulting selection is equivalent to plain set of indexes "index i1 i2 i3..."
    /// Useful to avoid recomputing selection on frame change when tracking given set of atoms
    void flatten();

    /// Returns a string formatted as Gromacs ndx file containing the current selection with given name.
    /// \warning Indexex in Gromacs ndx are starting from 1! Thus one is added to all pteros indexes!
    std::string to_gromacs_ndx(std::string name);

    friend void copy_coord(const Selection& from, int from_fr, Selection& to, int to_fr);        

    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Splitting selection into parts
    /// @{

    /// Split current selection into several selections according to
    /// the interatomic distances. Each resulting selection is a group
    /// of atoms connected by distances less than d.
    /// if d==0 the bonds from topology are used instead and periodic argument is ignored
    void split_by_connectivity(float d, std::vector<Selection>& res, bool periodic=true);    

    /// Split selection by residue index
    void split_by_residue(std::vector<Selection>& res);

    /// Split selection by molecule definition from topology
    void split_by_molecule(std::vector<Selection>& res);

    /// Split selection by chain
    void split_by_chain(std::vector<Selection>& chains);

    /// Split selection into contiguous ranges of indexes
    void split_by_contiguous_index(std::vector<Selection>& parts);

    /// Split selection into contiguous ranges of resindexes
    void split_by_contiguous_residue(std::vector<Selection>& parts);

    /// Split selection by the values returned by custom callback function.
    /// Callback is called for each atom in selection (local index is passed as its second argument)    
    /// Callback can return any type usable as a key of std::map.
    /// Selection is split into parts according to returned value.
    template<class F>
    void split(std::vector<Selection>& parts, F callback){
        using namespace std;
        parts.clear();
        using Ret = decltype(callback(*this,0));
        // Map
        map<Ret,vector<int> > m;
        for(int i=0; i<size(); ++i){
            m[callback(*this,i)].push_back(Index(i));
        }
        // Create selections
        typename map<Ret,vector<int> >::iterator it;
        for(it=m.begin();it!=m.end();it++){
            parts.push_back(Selection(*system));
            parts.back().modify( it->second.begin(), it->second.end() );
        }
    }


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
    *   it is possible to assign value to them. For example X(10) = 3.14
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

    /// Returns pointer to the coordinates of given atom for current frame.
    /// Used internally in Grid_searcher.
    inline Eigen::Vector3f* XYZ_ptr(int ind) const {
        return &(system->traj[frame].coord[index[ind]]);
    }

    /// Extracts X,Y and Z for given frame fr
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

    /// Extracts atomic number
    inline int& Element_number(int ind){
        return system->atoms[index[ind]].element_number;
    }

    inline const int& Element_number(int ind) const {
        return system->atoms[index[ind]].element_number;
    }

    /// Computes VDW radius. Read only.
    /// If atomic number is set uses whole periodic table from VMD
    /// If no atomic number uses rough guess from Gromacs vdwradii.dat
    float VDW(int ind) const;

    /// Extracts element name based on element_number. Read only.
    std::string Element_name(int ind) const;

    /** Returns periodic box of the frame pointed by selection
        The same as:
        \code
        sel.get_system()->Box(sel.get_frame());
        \endcode
        This is a convenience method. The same box is returnedby all selection
        which point to the same frame.
    */
    inline Periodic_box& Box() {
        return system->traj[frame].box;
    }

    inline const Periodic_box& Box() const {
        return system->traj[frame].box;
    }

    /** Returns time stamp of the frame pointed by selection
        The same as:
        \code
        sel.get_system()->Time(sel.get_frame());
        \endcode
        This is a convenience method. The same time is returnedby all selection
        which point to the same frame.
    */
    inline float& Time() {
        return system->traj[frame].time;
    }

    inline const float& Time() const {
        return system->traj[frame].time;
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
    std::unique_ptr<Selection_parser> parser;
    void allocate_parser();
    void sort_and_remove_duplicates();
};

//==============================================================================

/// Auxilary type used to incapsulate the atom and its current coordinates
/// Used internally in Selection::operator[] and in iterator access to Selection.
/// Objects of this class should not be created by the user in normal situation.
class Atom_proxy {
    friend class Selection::iterator;
    friend class Selection;
public:    
    const Selection& selection(){ return *sel; }
    const System& system(){ return *(sel->get_system()); }

    /// @name Inline accessors. Const and non-const versions.
    /// @{
    inline int& resid(){ return sel->Resid(ind); }
    inline const int& resid() const { return sel->Resid(ind); }

    inline const int& index() const { return sel->Index(ind); }

    inline std::string& name(){ return sel->Name(ind); }
    inline const std::string& name() const { return sel->Name(ind); }

    inline char& chain(){ return sel->Chain(ind); }
    inline const char& chain() const { return sel->Chain(ind); }

    inline std::string& resname(){ return sel->Resname(ind); }
    inline const std::string& resname() const { return sel->Resname(ind); }

    inline std::string& tag(){ return sel->Tag(ind); }
    inline const std::string& tag() const { return sel->Tag(ind); }

    inline float& occupancy(){ return sel->Occupancy(ind); }
    inline const float& occupancy() const { return sel->Occupancy(ind); }

    inline float& beta(){ return sel->Beta(ind); }
    inline const float& beta() const { return sel->Beta(ind); }

    inline int& resindex(){ return sel->Resindex(ind); }
    inline const int& resindex() const { return sel->Resindex(ind); }

    inline float& mass(){ return sel->Mass(ind); }
    inline const float& mass() const { return sel->Mass(ind); }

    inline float& charge(){ return sel->Charge(ind); }
    inline const float& charge() const { return sel->Charge(ind); }

    inline int& type(){ return sel->Type(ind); }
    inline const int& type() const { return sel->Type(ind); }

    inline std::string& type_name(){ return sel->Type_name(ind); }
    inline const std::string& type_name() const { return sel->Type_name(ind); }

    inline float& x(){ return sel->X(ind); }
    inline const float& x() const { return sel->X(ind); }

    inline float& y(){ return sel->Y(ind); }
    inline const float& y() const { return sel->Y(ind); }

    inline float& z(){ return sel->Z(ind); }
    inline const float& z() const { return sel->Z(ind); }

    inline Eigen::Vector3f& xyz(){ return sel->XYZ(ind); }
    inline const Eigen::Vector3f& xyz() const { return sel->XYZ(ind); }

    inline float& x(int fr){ return sel->X(ind,fr); }
    inline const float& x(int fr) const { return sel->X(ind,fr); }

    inline float& y(int fr){ return sel->Y(ind,fr); }
    inline const float& y(int fr) const { return sel->Y(ind,fr); }

    inline float& z(int fr){ return sel->Z(ind,fr); }
    inline const float& z(int fr) const { return sel->Z(ind,fr); }

    inline Eigen::Vector3f& xyz(int fr){ return sel->XYZ(ind,fr); }
    inline const Eigen::Vector3f& xyz(int fr) const { return sel->XYZ(ind,fr); }

    inline int& element_number(){ return sel->Element_number(ind); }
    inline const int& element_number() const { return sel->Element_number(ind); }

    std::string element_name() const { return sel->Element_name(ind); }
    float vdw() const { return sel->VDW(ind); }

    inline Atom& atom(){ return sel->Atom_data(ind); }
    inline const Atom& atom() const { return sel->Atom_data(ind); }

    /// @}

    /// Equality operator
    bool operator==(const Atom_proxy& other) const {
        return (sel==other.sel && ind==other.ind);
    }

    /// Inequality operator
    bool operator!=(const Atom_proxy &other) const {
        return !(*this == other);
    }    

private:
    Atom_proxy(){}
    Atom_proxy(Selection* s, int i): sel(s), ind(i) {}    

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
    iterator operator++(int junk) { iterator tmp = *this; ++proxy.ind; return tmp; }
    iterator& operator++() { ++proxy.ind; return *this; }
    Atom_proxy& operator*() { return proxy; }
    Atom_proxy* operator->() { return &proxy; }
    bool operator==(const iterator& rhs) { return proxy == rhs.proxy; }
    bool operator!=(const iterator& rhs) { return proxy != rhs.proxy; }
private:
    Atom_proxy proxy;
};

} // namespace pteros
#endif /* SELECTION_H */
