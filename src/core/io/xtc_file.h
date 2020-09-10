/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2020, Semen Yesylevskyy
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

#ifdef USE_GROMACS
#include "gromacs/fileio/xtcio.h"
#else
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#endif


namespace pteros {


class XTC_file: public Mol_file {
public:
    XTC_file(std::string& fname): Mol_file(fname), handle(nullptr) {}
    virtual void open(char open_mode);
    virtual ~XTC_file();

    virtual Mol_file_content get_content_type() const {        
        return Mol_file_content().traj(true).rand(true);
    }

protected:        

    virtual void do_write(const Selection &sel, const Mol_file_content& what) override;
    virtual bool do_read(System *sys, Frame *frame, const Mol_file_content& what) override ;
#ifdef USE_GROMACS
    virtual void seek_frame(int fr) override;
    virtual void seek_time(float t) override;
    virtual void tell_current_frame_and_time(int& step, float& t) override;
    virtual void tell_last_frame_and_time(int& step, float& t) override;
#endif

private:
#ifdef USE_GROMACS
    // for gmxlib
    t_fileio* handle;
    matrix box;
    int64_t step;
    int steps_per_frame;
    int64_t num_frames;
    float dt, max_t;
#else
    // for xdrfile
    XDRFILE* handle;
    matrix box;
    int step;
#endif
};

}



