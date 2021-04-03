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

#include "pteros/core/file_handler.h"
#include "xdrfile.h"
#include "xdrfile_trr.h"

namespace pteros {


class TrrFile: public FileHandler {
public:
    TrrFile(std::string& fname): FileHandler(fname), handle(nullptr) {}
    virtual void open(char open_mode);
    virtual ~TrrFile();

    virtual FileContent get_content_type() const {
        return FileContent()
                .traj(true);
    }

protected:

    virtual void do_write(const Selection &sel, const FileContent& what);
    virtual bool do_read(System *sys, Frame *frame, const FileContent& what);

private:
    // for xdrfile
    XDRFILE* handle;
    matrix box;
    int step;
};

}





