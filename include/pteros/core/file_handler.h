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
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include <bitset>

namespace pteros {

class FileContent {
public:
    FileContent(){
        flags.reset();
    }

    // The list of atoms and their properties
    bool atoms() const { return flags[0]; }
    FileContent atoms(bool val){ flags[0] = val; return *this;}

    // Single set of coordinates
    bool coord() const { return flags[1]; }
    FileContent coord(bool val){ flags[1] = val; return *this;}

    // Many frames
    bool traj() const { return flags[2]; }
    FileContent traj(bool val){ flags[2] = val; return *this;}

    // Molecular topology
    bool top() const { return flags[3]; }
    FileContent top(bool val){ flags[3] = val; return *this;}

    // Random access trajectory
    bool rand() const { return flags[4]; }
    FileContent rand(bool val){ flags[4] = val; return *this;}

private:
    std::bitset<5> flags;
};

//------------------------------------------------------------------------------

class FileHandler;
using FileHandler_ptr = std::unique_ptr<FileHandler>;

/// Generic API for reading and writing any molecule file formats
class FileHandler {
public:
    /// Recognizes file extension and returns a file handler
    static FileHandler_ptr recognize(std::string fname);

    /** Recognize file extension, open file for reading or writing and return a file handler.
     This function is aquivalent to:
     \code
     auto f = Mol_file::recognize(fname);
     f.open(mode);
     \endcode
    */
    static FileHandler_ptr open(std::string fname, char open_mode);

    /// Opens a file with given access mode. Need to be defined by derived classes.
    virtual void open(char open_mode) = 0;
    virtual void close();

    virtual ~FileHandler();

    /// Reads data, which are specified by what.
    /// Pointers to System and Frame could be nullptr if not used
    /// Returns true if read operation is succesfull and false if not.    
    bool read(System* sys, Frame* frame, const FileContent& what);

    /// Write data from selection specidied by what.
    void write(const Selection& sel, const FileContent& what);

    /// Reports content of this file type
    virtual FileContent get_content_type() const = 0;

protected:    
    FileHandler(std::string& file_name);

    // Stores file name
    std::string fname;
    // Number of atoms
    int natoms;    

    // Method to sanity check parameters send to read and write
    void sanity_check_read(System* sys, Frame* frame, const FileContent &what) const;
    void sanity_check_write(const Selection& sel, const FileContent& what) const;

    /// User-overriden method for reading
    virtual bool do_read(System* sys, Frame* frame, const FileContent& what) = 0;

    /// User-overriden method for writing
    virtual void do_write(const Selection& sel, const FileContent& what) = 0;
};

//------------------------------------------------------------------------------

class FileHandlerRandomAccess: public FileHandler {
protected:
    FileHandlerRandomAccess(std::string& file_name): FileHandler(file_name) {}

public:
    // Seek frame
    virtual void seek_frame(int fr) = 0;

    // Seek time
    virtual void seek_time(float t) = 0;

    // Report current position in trajectory
    virtual void tell_current_frame_and_time(int& step, float& t) = 0;

    // Report last position in trajectory
    virtual void tell_last_frame_and_time(int& step, float& t) = 0;
};

//------------------------------------------------------------------------------

} // namespace pteros


