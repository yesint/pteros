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

namespace pteros {

class XtcFile: public FileHandlerRandomAccess {
public:    
    XtcFile(const std::string& fname): FileHandlerRandomAccess(fname), content(FileContent().traj(true).rand(true)) {}
    virtual void open(char open_mode) override;
    virtual void close() override;
    virtual ~XtcFile();

    virtual FileContent get_content_type() const {
        return content;
    }

protected:

    virtual void do_write(const Selection &sel, const FileContent& what) override;
    virtual bool do_read(System *sys, Frame *frame, const FileContent& what) override ;

    virtual void seek_frame(int fr) override;
    virtual void seek_time(float t) override;
    virtual void tell_current_frame_and_time(int& step, float& t) override;
    virtual void tell_last_frame_and_time(int& step, float& t) override;

private:
    FileContent content;

    struct XTC_internals;
    XTC_internals* xtc;
};

}





