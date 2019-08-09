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


#pragma once

#include "pteros/core/mol_file.h"
#include "gromacs/fileio/xtcio.h"

namespace pteros {


class XTC_file: public Mol_file {
public:
    XTC_file(std::string& fname): Mol_file(fname) {}
    virtual void open(char open_mode);
    virtual ~XTC_file();

    virtual Mol_file_content get_content_type() const {        
        return Mol_file_content().traj(true);
    }

protected:        

    virtual void do_write(const Selection &sel, const Mol_file_content& what);

    virtual bool do_read(System *sys, Frame *frame, const Mol_file_content& what);

private:
    t_fileio* handle;
    matrix box;
    int64_t step;
    bool first;
};

}


