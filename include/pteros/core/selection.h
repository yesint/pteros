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
//#include <memory>
#include <map>
#include <vector>
#include <functional>
#include <Eigen/Core>
//#include <Eigen/Geometry>
#include "pteros/core/atom_handler.h"
#include "pteros/core/system.h"
#include "pteros/core/typedefs.h"

namespace pteros {

// Forward declaration of friend classes
class System;
class SelectionParser;


/** @brief Selection class.
*
*   Selections are key objects in Pteros. Technically speaking the selection
*   is just an array, which contains indexes of selected atoms in particular system.
*   Selection does not hold the copies of the atoms or their coordinates, it
*   just points to them serving like an alias for certain subset of atoms.
*   Thus selections may overlap arbitrarily.
*   Selections are used to perform various operations on the group of selected atoms.
*   The changes become immediately visible to all other selections, which point to
*   some of changed atoms.
*   Each selection is bound to particular System. There are neither 'parentless' selection nor the
*   selections, which combine the atoms from different systems.
*   Selections are created using the syntax, which is very similar to those used in VMD but with
*   many useful extensions.
*/
class Selection {
  friend class System;
  friend class SelectionParser;
  friend class Grid_searcher;

  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Constructors and operators
    /// @{

    /// Default constructor for empty selection.
    Selection();

    /** Associates selection with the system @param sys,
    *   Selection content should be set later by calling modify().
    */
    explicit Selection(const System& sys);

    /** Textual selection constructor.
    *   @param sys System pointed by this selection
    *   @param str Selection string    
    */
    Selection(const System& sys, std::string str, int fr = 0);

    /** Constructor, which creates selection from the interval of indexes.
        It is much faster then parsing selection corresponding string, but is limited
        to contigous interval of indexes.
        @param sys System pointed by this selection
        @param ind1 First index in interval
        @param ind2 Last index in interval (inclusive!)        
     */
    Selection(const System& sys, int ind1, int ind2);

    /// Constructor from the vector of indexes
    /// Vector may be in any order and may contain duplicates.
    Selection(const System& sys, const std::vector<int>& ind);

    /// Constructor from the pair of iterators to int sequence
    /// Sequence may be in any order and may contain duplicates.
    Selection(const System& sys,
              std::vector<int>::iterator it1,
              std::vector<int>::iterator it2);

    /** Constructor which takes user-defined callback function.
      Callback takes the system as first argument, target frame number as the second
      and the vector to be filled by selected atom indexes.
      Vector may be filled in any order and may contain duplicates.
      \warning
      Resulting selection is neither coordinate-dependent nor text-based.
      It won't recompute itself on the frame change even if it involves coordinates.
    */
    Selection(const System& sys,
              const std::function<void(const System&,int,std::vector<int>&)>& callback,
              int fr = 0);

    /// Copy constructor
    Selection(const Selection& other);

    /// Destructor
    virtual ~Selection();

    /// Assignment operator
    Selection& operator=(const Selection& other);

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
     sel[0].name() = "A"; // OK to write fields of proxy object
     \endcode

     Main goal of [] operator is usage in range-based for loops:
     \code
     Selection sel(sys,"name CA");
     for(auto a: sel){
        cout << a.name() << " " << a.x() << endl;
     }
     \endcode

     Otherwise it is slower than conventional syntax like sel.name(i)
     */
    AtomHandler operator[](int ind);

    /** Indexing operator which sets both index and frame.
     * Takes a pair of values.
      \code
      Selection sel(sys,"name CA");
      int at = 123;
      // Assign x coorditate of atom at for frame 0 from frame 1.
      sel[{at,0}].x() = sel[{at,1}].x();
      \endcode

    */
    AtomHandler operator[](const std::pair<int,int>& ind_fr);

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

    /// Append several selections to this one
    void append(const std::vector<Selection>& sel_vec);

    /// Remove all atoms of sel from current selection
    void remove(const Selection& sel);

    /// Remove given absolute index from current selection
    void remove(int ind);

    /// Inverts selection in place by selecting those atoms which were not selected
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

    /// Convenience method, which combines set_system and modify(str)
    void modify(const System& sys, std::string str, int fr = 0);

    /// Convenience method, which combines set_system and modify(int,int)
    void modify(const System& sys, int ind1, int ind2);

    /// Convenience method, which combines set_system and modify(vector<int>)
    void modify(const System& sys, const std::vector<int>& ind);

    /// Convenience method, which combines set_system and modify(iter,iter)
    void modify(const System& sys,
                std::vector<int>::iterator it1,
                std::vector<int>::iterator it2);

    /// Convenience method, which combines set_system and modify(callback)
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
    std::vector<int>::const_iterator index_begin() const { return _index.begin(); }

    /// Get const iterator for the end of index
    std::vector<int>::const_iterator index_end() const { return _index.end(); }

    /// Get back_insert_iterator for index
    std::back_insert_iterator<std::vector<int> > index_back_inserter() {
        return std::back_inserter(_index);
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
    std::vector<int> get_index() const { return _index; }

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


    /// Get charge
    std::vector<float> get_charge() const;

    /// Set charge in selection to the values from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_charge(const std::vector<float>& data);

    /// Sets charge of all selected atoms to the same given value.
    void set_charge(float data);

    /// Computes total charge of selection
    float get_total_charge() const;

    /// Get tags
    std::vector<std::string> get_tag(bool unique=false) const;

    /// Set tags in selection to the values from supplied vector.
    /// \note Vector size must be the save as the size of selection.
    void set_tag(const std::vector<std::string>& data);

    /// Set tags of all selected atoms to the same given value.
    void set_tag(std::string data);


    /// Get velocities of all atoms in this selection for the current frame
    /// By default columnt of the matrix contain atom coordinates.
    /// Optional flag turns this to row-major format (mainly used by python binding)
    Eigen::MatrixXf get_vel(bool make_row_major_matrix = false) const;

    /// Set velocities of this selection for current frame
    void set_vel(MatrixXf_const_ref data);

    /// Get forces of all atoms in this selection for the current frame
    /// By default columnt of the matrix contain atom coordinates.
    /// Optional flag turns this to row-major format (mainly used by python binding)
    Eigen::MatrixXf get_force(bool make_row_major_matrix = false) const;

    /// Set forces of this selection for current frame
    void set_force(MatrixXf_const_ref coord);

    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Computing properties of selection
    /// @{

    /// Checks if selection is larger than 1/2 of the box size in
    /// any dimension. Returns true if so.
    bool is_large();

    /** Get the center of selection.
     @param mass_weighted Use mass-weighting
     @param periodic Account for periodic boundary conditions.

     \warning
     If the size of selection is larger than 1/2 of the box size in
     any dimension you will get incorrect results if periodic is set to true.
     This is not checked automatically!
     In this case use one of unwrapping options first.
    */
    Eigen::Vector3f center(bool mass_weighted = false,
                           Array3i_const_ref pbc = noPBC,
                           int pbc_atom = -1) const;

    // Get the center of selection with user-supplied weights
    Eigen::Vector3f center(const std::vector<float>& weights,
                           Array3i_const_ref pbc = noPBC,
                           int pbc_atom = -1) const;

    /// Get minimal and maximal coordinates in selection
    void minmax(Vector3f_ref min, Vector3f_ref max) const;

    /// Get the SASA using powersasa algorithm. Returns area and computes volume and per-atom values if asked
    float powersasa(float probe_r = 0.14,
               std::vector<float>* area_per_atom = nullptr,
               float* total_volume = nullptr,
               std::vector<float>* volume_per_atom = nullptr) const;

    /// Get SASA using Shrake and Rupley algorithm (slower than powersasa and can't compute volumes)
    float sasa(float probe_r = 0.14, std::vector<float>* area_per_atom = nullptr, int n_sphere_points = 960) const;

    /// Computes average structure over the range of frames
    Eigen::MatrixXf average_structure(int b=0, int e=-1, bool make_row_major_matrix = false) const;

    /** Extracts X,Y,Z for given atom index for specified range of frames
    *   (gets trajectory of given atom).
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
    void inertia(Vector3f_ref moments, Matrix3f_ref axes,
                 Array3i_const_ref pbc = noPBC,
                 int pbc_atom = -1) const;

    /** Computes radius of gyration for selection
     \warning
     If the size of selection is larger than 1/2 of the box size in
     any dimension you will get incorrect results if periodic is set to true.
     This is not checked automatically!
     In this case use one of unwrapping options first.
    */
    float gyration(Array3i_const_ref pbc = noPBC, int pbc_atom = -1) const;

    /* Computes dipole moment of selection in Debye
    If the size of selection is larger than 1/2 of the box size in
    any dimension you will get incorrect results if periodic is set to true.
    This is not checked automatically!
    In this case use one of unwrapping options first.
    */
    Eigen::Vector3f dipole(bool as_charged=false, Array3i_const_ref pbc = fullPBC, int pbc_atom = -1) const;

    /// Get distance between two atoms (periodic in given dimensions if needed).
    /// \note
    /// This function takes selection indexes, not absolute indexes.
    float distance(int i, int j, Array3i_const_ref pbc = fullPBC) const;

    /// Get angle in degrees between three atoms (periodic in given dimensions if needed).
    /// \note
    /// This function takes selection indexes, not absolute indexes.
    float angle(int i, int j, int k, Array3i_const_ref pbc = fullPBC) const;

    /// Get dihedral angle in degrees between three atoms (periodic in given dimensions if needed).
    /// \note
    /// This function takes selection indexes, not absolute indexes.
    float dihedral(int i, int j, int k, int l, Array3i_const_ref pbc = fullPBC) const;
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Geometry transformation functions
    /// @{

    /// Translate selection by given vector
    void translate(Vector3f_const_ref v);

    /// Translate center of masses of selection to given point
    void translate_to(Vector3f_const_ref p,
                      bool mass_weighted = false,
                      Array3i_const_ref pbc = noPBC,
                      int pbc_atom = -1);

    /// Rotate selection around given axis relative to given pivot
    /// @param pivot Rotation around this pivot
    /// @param axis Rotate around this vector
    /// @param angle Rotation angle in radians
    void rotate(Vector3f_const_ref pivot, Vector3f_const_ref axis, float angle);

    /// Wraps whole selection to the periodic box
    void wrap(Array3i_const_ref pbc = fullPBC);

    /** Unwraps selection to make it whole if possible (without jumps over periodic box boundary).
      The periodic center of mass is used as an anchor point if pbc_atom<0.
      \warning
      If the size of selection is larger than 1/2 of the box size in
      any dimension unwrap() will not work as expected and will not make selection "compact"!
      This is \e not checked automatically.
    */
    void unwrap(Array3i_const_ref pbc = fullPBC, int pbc_atom = -1);

    /** Unwraps selection to make it whole (without jumps over periodic box boundary).
     * based on preserving all bonds.
     * This method works reliably in any case, but is much slower than unwrap()
     * @param d Maximal bond length. If 0 bonds from topology are used (if present).
     * @param pbc_atom Local index of the reference atom, which doesn't move.
     * @return Number of disconnected pieces after unwrapping. 1 means solid selection.
     */
    int unwrap_bonds(float d, Array3i_const_ref pbc = fullPBC, int pbc_atom = -1);

    /** Get transform for orienting selection by principal axes.
     * \warning
     * If the size of selection is larger than 1/2 of the box size in
     * any dimension you will get funny results if @param is_periodic is set to true.
     */    
    Eigen::Affine3f principal_transform(Array3i_const_ref pbc = noPBC, int pbc_atom = -1) const;

    /** Orient molecule by its principal axes.
     * The same as
     * \code
     * Eigen::Affine3f tr = sel.principal_transform();
     * sel.apply_transform(tr);
     * \endcode
     */
    void principal_orient(Array3i_const_ref pbc = noPBC, int pbc_atom = -1);

    /// Returns total number of residues in selection. Partial residues are also included.
    int num_residues() const;
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

    /// Apply transformation
    void apply_transform(const Eigen::Affine3f& t);
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Energy functions
    /// @{

    /// Self-energy of selection computed within given interaction cut-off.
    /// If cutoff is 0 the cutoff from topology is used.
    Eigen::Vector2f non_bond_energy(float cutoff=0, bool pbc = true) const;


    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name File IO
    /// @{

    /** Write structure or trajectory for selection. File type is determined by extension.
    *   Frames from b to e are written.
    *   If @param b is not set or -1 it means current frame
    *   If @param e is not set or -1 it means the last frame
    */
    void write(std::string fname, int b=-1,int e=-1) const;

    void write(const std::unique_ptr<FileHandler>& handler, FileContent what,int b=-1,int e=-1) const;
    /// @}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Utility functions
    /// @{

    /// Get the size of selection
    int size() const {return _index.size();}

    /// Returns true if selection was created from text string and false if it was
    /// constructed 'by hand' by appending indexes or other selections
    bool text_based() const;

    /// Returns true if selection is coordinate-dependent and is able to recompute
    /// itself on the change of frame
    bool coord_dependent() const;

    /// "Flattens" selection by removing coordinate dependence and making it not text-based.
    /// Resulting selection is equivalent to plain set of indexes "index i1 i2 i3..."
    /// Useful to avoid recomputing selection on frame change when tracking given set of atoms
    void flatten();

    /// Returns a string formatted as Gromacs ndx file containing the current selection with given name.
    /// \warning Indexes in Gromacs ndx are starting from 1! Thus one is added to all pteros indexes!
    std::string to_gromacs_ndx(std::string name) const;

    /// Copy coordinates from one selection to the other using given frames.
    /// \warning Size of selections should be the same.
    friend void copy_coord(const Selection& from, int from_fr, Selection& to, int to_fr);        

    /// Copy coordinates from one selection to the other using their current frames.
    /// \warning Size of selections should be the same.
    friend void copy_coord(const Selection& from, Selection& to);

    /// Finds local index in selection of provided global index
    /// Returns -1 if not found
    int find_index(int global_index) const;

    /// Get all bonds within this selection returned as local selection indexes in form 1->[2,3,..].
    /// If d>0 it is used as cutoff
    /// if d==0 the bonds from topology are used and periodicity is ignored    
    std::vector<std::vector<int>> get_internal_bonds(float d, bool periodic=true) const;
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
    /// Signature of callback is:
    /// T cb(const Selection&,int)
    /// Callback is called for each atom in selection (local index is passed as its second argument)    
    /// Callback can return any type usable as a key of std::map.
    /// Selection is split into parts according to returned values - all identical values
    /// go to the same part.
    template<class F>
    void split(std::vector<Selection>& parts, F callback){
        using namespace std;
        parts.clear();
        using Ret = decltype(callback(*this,0));
        map<Ret,vector<int> > m;
        for(int i=0; i<size(); ++i){
            m[callback(*this,i)].push_back(index(i));
        }
        // Create selections
        typename map<Ret,vector<int> >::iterator it;
        for(it=m.begin();it!=m.end();it++){
            parts.emplace_back(*system, it->second.begin(), it->second.end());
        }
    }


    /// Selects each residue, which is referenced by selection.
    /// All selections for residues are placed into supplied vector.
    /// Handles multiple chains correctly.    
    void each_residue(std::vector<Selection>& sel) const;
    /// @}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// @name Secondary structure methods
    /// These methods take protein residues referenced by selection and compute DSSP on them.
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
    /** @name Inline accessors
    *   Used to access the properties of
    *   particular atom in selection. The ind passed to these functions is
    *   is the selection index, not the global system %index. I.e. passing 10 will
    *   extract the property of %atom with global %index index[10]
    *   <br>All these function except index could be used as lvalue, which means that
    *   it is possible to assign value to them. For example x(10) = 3.14
    */
    /// @{

    /// Extracts X for current frame
    inline float& x(int ind)       { return system->traj[frame].coord[_index[ind]](0);  }
    inline float  x(int ind) const { return system->traj[frame].coord[_index[ind]](0); }
    /// Extracts X for given frame frame fr
    inline float& x(int ind, int fr)       { return system->traj[fr].coord[_index[ind]](0); }
    inline float  x(int ind, int fr) const { return system->traj[fr].coord[_index[ind]](0); }

    /// Extracts Y for current frame
    inline float& y(int ind)       { return system->traj[frame].coord[_index[ind]](1); }
    inline float  y(int ind) const { return system->traj[frame].coord[_index[ind]](1); }
    /// Extracts Y for given frame frame fr
    inline float& y(int ind, int fr)       { return system->traj[fr].coord[_index[ind]](1); }
    inline float  y(int ind, int fr) const { return system->traj[fr].coord[_index[ind]](1); }

    /// Extracts Z for current frame
    inline float& z(int ind)       { return system->traj[frame].coord[_index[ind]](2); }
    inline float  z(int ind) const { return system->traj[frame].coord[_index[ind]](2); }
    /// Extracts Z for given frame frame fr
    inline float& z(int ind, int fr)       { return system->traj[fr].coord[_index[ind]](2); }
    inline float  z(int ind, int fr) const { return system->traj[fr].coord[_index[ind]](2); }

    /// Extracts X,Y and Z for current frame
    inline       Eigen::Vector3f& xyz(int ind)       { return system->traj[frame].coord[_index[ind]]; }
    inline const Eigen::Vector3f& xyz(int ind) const { return system->traj[frame].coord[_index[ind]]; }
    /// Extracts X,Y and Z for given frame fr
    inline       Eigen::Vector3f& xyz(int ind, int fr)       { return system->traj[fr].coord[_index[ind]]; }
    inline const Eigen::Vector3f& xyz(int ind, int fr) const { return system->traj[fr].coord[_index[ind]]; }

    /// Returns pointer to the coordinates of given atom for current frame.
    /// Used internally in Grid_searcher.
    inline Eigen::Vector3f* xyz_ptr(int ind)         const { return &(system->traj[frame].coord[_index[ind]]); }
    inline Eigen::Vector3f* xyz_ptr(int ind, int fr) const { return &(system->traj[fr].coord[_index[ind]]); }

    /// Extracts velocity for current frame
    inline       Eigen::Vector3f& vel(int ind)       { return system->traj[frame].vel[_index[ind]]; }
    inline const Eigen::Vector3f& vel(int ind) const { return system->traj[frame].vel[_index[ind]]; }
    /// Extracts velocity for given frame fr
    inline       Eigen::Vector3f& vel(int ind, int fr)       { return system->traj[fr].vel[_index[ind]]; }
    inline const Eigen::Vector3f& vel(int ind, int fr) const { return system->traj[fr].vel[_index[ind]]; }

    /// Extracts force for current frame
    inline       Eigen::Vector3f& force(int ind)       { return system->traj[frame].force[_index[ind]]; }
    inline const Eigen::Vector3f& force(int ind) const { return system->traj[frame].force[_index[ind]]; }
    /// Extracts force for given frame fr
    inline       Eigen::Vector3f& force(int ind, int fr)       { return system->traj[fr].force[_index[ind]]; }
    inline const Eigen::Vector3f& force(int ind, int fr) const { return system->traj[fr].force[_index[ind]]; }

#define DEFINE_ACCESSOR(T,prop) \
    inline T& prop(int ind){ return system->atoms[_index[ind]].prop; } \
    inline const T& prop(int ind) const { return system->atoms[_index[ind]].prop; }

    /// Extracts type
    DEFINE_ACCESSOR(int,type)
    /// Extracts typename
    DEFINE_ACCESSOR(std::string,type_name)
    /// Extracts residue name
    DEFINE_ACCESSOR(std::string,resname)
    /// Extracts chain
    DEFINE_ACCESSOR(char,chain)
    /// Extracts atom name
    DEFINE_ACCESSOR(std::string,name)
    /// Extracts atom mass
    DEFINE_ACCESSOR(float,mass)
    /// Extracts atom charge
    DEFINE_ACCESSOR(float,charge)
    /// Extracts B-factor
    DEFINE_ACCESSOR(float,beta)
    /// Extracts occupancy field
    DEFINE_ACCESSOR(float,occupancy)
    /// Extracts residue number
    DEFINE_ACCESSOR(int,resid)
    /// Extracts tag
    DEFINE_ACCESSOR(std::string,tag)

    /// Extracts atom index in the system, which is pointed by selection
    inline int index(int ind) const { return _index[ind]; }

    /// Extracts whole atom
    inline Atom& atom(int ind){ return system->atoms[_index[ind]]; }
    inline const Atom& atom(int ind) const { return system->atoms[_index[ind]]; }

    /// Extracts resindex
    DEFINE_ACCESSOR(int,resindex)
    /// Extracts atomic number
    DEFINE_ACCESSOR(int,atomic_number)

    /// Computes VDW radius. Read only.
    /// If atomic number is set uses whole periodic table from VMD
    /// If no atomic number uses rough guess from Gromacs vdwradii.dat
    float vdw(int ind) const;

    /// Extracts element name based on element_number. Read only.
    std::string element_name(int ind) const;

    /** Returns periodic box of the frame pointed by selection
        The same as:
        \code
        sel.get_system()->box(sel.get_frame());
        \endcode
        This is a convenience method. The same box is returned by all selection
        which point to the same frame.
    */
    inline PeriodicBox& box() { return system->traj[frame].box; }
    inline const PeriodicBox& box() const { return system->traj[frame].box; }

    /** Returns time stamp of the frame pointed by selection
        The same as:
        \code
        sel.get_system()->Time(sel.get_frame());
        \endcode
        This is a convenience method. The same time is returned by all selection
        which point to the same frame.
    */
    inline float& time() { return system->traj[frame].time; }
    inline const float& time() const { return system->traj[frame].time; }

    /// @}

protected:
    // Raw text of selection
    std::string sel_text;    
    // Indexes of atoms in selection
    std::vector<int> _index;
    // Pointer to target system
    System* system;

    // Stores current frame
    int frame;

    // Holds an instance of selection parser
    std::unique_ptr<SelectionParser> parser;
    void allocate_parser();
    void sort_and_remove_duplicates();    
    void process_pbc_atom(int& a) const;
    void get_local_bonds_from_topology(std::vector<std::vector<int>>& con) const;
};

//-----------------------------------------------------------------------
/// Random-access forward iterator for Selection
class Selection::iterator {
public:
    using value_type = AtomHandler;
    using difference_type = size_t;
    using pointer = AtomHandler*;
    using reference = AtomHandler&;
    using iterator_category = std::random_access_iterator_tag;

    iterator(const Selection& sel, int i): ind(i) {
        proxy.set(sel,ind);
    }
    iterator operator++(int junk) { iterator tmp = *this; ++ind; return tmp; }
    iterator& operator++() { ++ind; proxy.advance(); return *this; }
    iterator& operator+(int i) {ind+=i; proxy.advance(i); return *this;}
    iterator& operator-(int i) {ind-=i; proxy.advance(-i); return *this;}
    reference operator*() { return proxy; }
    pointer   operator->() { return &proxy; }
    bool operator==(const iterator& rhs) { return ind == rhs.ind; }
    bool operator!=(const iterator& rhs) { return !(*this==rhs); }

private:
    int ind;
    AtomHandler proxy;
};

/// Checks if several selections overlap
bool check_selection_overlap(const std::vector<Selection> &sel_vec);


/// Non-bond energy between two selections computed within given interaction cut-off.
/// If cutoff is 0 the cutoff from topology is used.
/// fr = -1 computes for current frame of selection 1.
Eigen::Vector2f non_bond_energy(const Selection& sel1,
                                const Selection& sel2,
                                float cutoff = 0,
                                int fr = -1,
                                bool pbc = true);

} // namespace pteros





