/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2023, Semen Yesylevskyy
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
#include <cstdio>
#include "pteros/core/file_handler.h"

namespace pteros {

/// Reader for GRO files. It doesn't use VMD plugins because it doesn't support writing
class GroFile: public FileHandler {
public:
    // High-level API        
    GroFile(const std::string& fname): FileHandler(fname) {}
    virtual void open(char open_mode) override;
    virtual void close() override;

    virtual FileContent get_content_type() const override {        
        return FileContent()
                .atoms(true)
                .coord(true);
    }

protected:

    std::FILE* file_handle;

    virtual bool do_read(System *sys, Frame *frame, const FileContent& what) override;
    virtual void do_write(const Selection &sel, const FileContent& what) override;
};

}




