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


#ifndef MOL_FILE_H
#define MOL_FILE_H

#include <string>
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include <bitset>

namespace pteros {

class Mol_file_content {
public:
    Mol_file_content(){
        flags.reset();
    }

    // The list of atoms and their properties
    bool atoms() const { return flags[0]; }
    Mol_file_content atoms(bool val){ flags[0] = val; return *this;}

    // Single set of coordinates
    bool coord() const { return flags[1]; }
    Mol_file_content coord(bool val){ flags[1] = val; return *this;}

    // Many frames
    bool traj() const { return flags[2]; }
    Mol_file_content traj(bool val){ flags[2] = val; return *this;}

    // Molecular topology
    bool top() const { return flags[3]; }
    Mol_file_content top(bool val){ flags[3] = val; return *this;}

    // Random access trajectory (for the future)
    bool rand() const { return flags[4]; }
    Mol_file_content rand(bool val){ flags[4] = val; return *this;}

private:
    std::bitset<5> flags;
};

/// Generic API for reading and writing any molecule file formats
class Mol_file {
public:
    /// Recognizes file extension and returns a file handler
    static std::unique_ptr<Mol_file> recognize(std::string fname);

    /** Recognize file extension, open file for reading or writing and return a file handler.
     This function is aquivalent to:
     \code
     auto f = Mol_file::recognize(fname);
     f.open(mode);
     \endcode
    */
    static std::unique_ptr<Mol_file> open(std::string fname, char open_mode);

    /// Opens a file with given access mode. Need to be defined by derived classes.
    virtual void open(char open_mode) = 0;

    virtual ~Mol_file();

    /// Reads data, which are specified by what.
    /// Pointers to System and Frame could be nullptr if not used
    /// Returns true if read operation is succesfull and false if not.    
    bool read(System* sys, Frame* frame, const Mol_file_content& what);

    /// Write data from selection specidied by what.
    void write(const Selection& sel, const Mol_file_content& what);

    /// Reports content of this file type
    virtual Mol_file_content get_content_type() const = 0;

protected:    
    Mol_file(std::string& file_name);

    // Stores file name
    std::string fname;
    // Number of atoms
    int natoms;    
    // Functions called to update System on file reading
    // Mol_file is a friend of System and can access it's internals
    // but derived *_file classes are not friends and need to call these functions.
    void allocate_atoms_in_system(System& sys, int n);
    void set_atom_in_system(System& sys, int i, Atom& at);
    Atom& atom_in_system(System& sys, int i);
    void append_atom_in_system(System& sys, Atom& at);

    // Method to sanity check parameters send to read and write
    void sanity_check_read(System* sys, Frame* frame, const Mol_file_content &what) const;
    void sanity_check_write(const Selection& sel, const Mol_file_content& what) const;

    /// User-overriden method for reading
    virtual bool do_read(System* sys, Frame* frame, const Mol_file_content& what) = 0;

    /// User-overriden method for writing
    virtual void do_write(const Selection& sel, const Mol_file_content& what) = 0;
};

float get_mass_from_atom_name(std::string& name);

}
#endif /* MOL_FILE_H */

