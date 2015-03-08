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

#ifndef MOL_FILE_H
#define MOL_FILE_H

#include <string>
#include "pteros/core/system.h"
#include "pteros/core/selection.h"

namespace pteros {

/// What is stored in this file type?
struct Mol_file_content {
    bool structure; // The list of atoms and their properties
    bool coordinates; // Single frame
    bool trajectory; // Many frames
    bool topology; // Molecular topology

    Mol_file_content(): structure(false), coordinates(false),
                        trajectory(false), topology(false) {}

    Mol_file_content(bool str, bool coor, bool trj, bool top): structure(str),
        coordinates(coor), trajectory(trj), topology(top) {}
};

/// Generic API for reading and writing any molecule file formats
class Mol_file {
public:
    // Recognizes file extension and returns a handler object
    static std::unique_ptr<Mol_file> recognize(std::string fname);
    // Opens a file with given access mode. Need to be defined by derived classes.
    virtual void open(char open_mode) = 0;

    virtual ~Mol_file();

    /// Reads data, which are specified by what.
    /// Pointers to System and Frame could be NULL if not used
    bool read(System* sys, Frame* frame, const Mol_file_content& what);

    /// Write given data from selection specidied by what.
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
    Force_field& ff_in_system(System& sys);

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
