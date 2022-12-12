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

#include "xtc_file.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "gromacs_utils.h"

#include "xdr_utils.h"
#include "xdrfile.h"
#include "xdrfile_xtc.h"


using namespace std;
using namespace pteros;
using namespace Eigen;

// xdrfile internals
struct XtcFile::XTC_internals {
    XDRFILE* handle;
    matrix box;
    int step;
    int steps_per_frame;
    int64_t num_frames;
    float dt, max_t;

    XTC_internals(): handle(nullptr) {}
};


void XtcFile::open(char open_mode)
{
    xtc = new XTC_internals();

    bool bOk;
    xtc->handle = xdrfile_open(fname.c_str(),&open_mode);

    if(!xtc->handle) throw PterosError("Unable to open XTC file {}", fname);

    // Extract number of atoms
    int ok = xdr_xtc_get_natoms(xtc->handle,&natoms);
    if(!ok) throw PterosError("Can't read XTC number of atoms");

    if(open_mode=='r'){
        // XTC file contains step number in terms of simulation steps, not saved frames
        // So we have to extract conversion factor
        int next = xtc_get_next_frame_number(xtc->handle,natoms);
        int cur = xtc_get_current_frame_number(xtc->handle,natoms,&bOk);
        if(cur<0 || next<0 || !bOk) throw PterosError("Can't detect number of steps per frame");
        xtc->steps_per_frame = next-cur;

        // Get total number of frames in the trajectory
        xtc->num_frames = xdr_xtc_get_last_frame_number(xtc->handle,natoms,&bOk);

        if(!bOk) throw PterosError("Can't get number of frames");
        if(xtc->num_frames<0){
            LOG()->warn("Weird XTC file: negative number of frames returned ({})!",xtc->num_frames);
            LOG()->warn("Random access operations disabled on this trajectory.");
            content.rand(false);
        }
        xtc->num_frames /= xtc->steps_per_frame;

        // Get time step
        xtc->dt = xdr_xtc_estimate_dt(xtc->handle,natoms,&bOk);
        if(!bOk) throw PterosError("Can't get time step");

        xtc->max_t = xdr_xtc_get_last_frame_time(xtc->handle,natoms,&bOk);
        if(!bOk || xtc->max_t<0) throw PterosError("Can't get last frame time");

        LOG()->debug("There are {} frames, max_t= {}, dt={}",xtc->num_frames,xtc->max_t,xtc->dt);
    }

    if(!xtc->handle) throw PterosError("Unable to open XTC file {}", fname);

    // Prepare the box just in case
    init_gmx_box(xtc->box);

    // -1 for reading means initialization step
    xtc->step = (open_mode=='r') ? -1 : 0;
}

void XtcFile::close()
{
    if(xtc->handle) xdrfile_close(xtc->handle);
    if(xtc) delete xtc;
}

XtcFile::~XtcFile()
{
    close();
}

bool XtcFile::do_read(System *sys, Frame *frame, const FileContent &what){
    float prec;
    int ret;

    frame->coord.resize(natoms);
    ret = read_xtc(xtc->handle,natoms,&xtc->step,&frame->time,xtc->box, (rvec*)frame->coord.data(), &prec);
    if(ret == exdrENDOFFILE) return false; // End of file
    if(ret != exdrOK){
        LOG()->warn("XTC frame {} is corrupted!",xtc->step);
        return false;
    }

    if(xtc->step==0){
        LOG()->debug("Number of atoms: {}",natoms);
        LOG()->debug("XTC precision: {}",prec);
    }

    gmx_box_to_pteros(xtc->box,frame->box);
    return true;
}


void XtcFile::seek_frame(int fr)
{
    if(fr>=xtc->num_frames) throw PterosError("Can't seek to frame {}, there are {} frames in this file",fr,xtc->num_frames);
    int ret = xdr_xtc_seek_frame(fr*xtc->steps_per_frame,xtc->handle,natoms);
    if(ret<0) throw PterosError("Error seeking to frame {}",fr);
}

void XtcFile::seek_time(float t)
{
    if(t<0 || t>xtc->max_t) throw PterosError("Can't seek to time {}, last time is {}",t,xtc->max_t);
    // We assume equally spaced frames in the trajectory. It's much faster
    int ret = xdr_xtc_seek_frame(ceil(t/xtc->dt)*xtc->steps_per_frame,xtc->handle,natoms);
    //int ret = xdr_xtc_seek_time(t,handle,natoms,false);
    if(ret<0) throw PterosError("Can't seek to time {}",t);
}

void XtcFile::tell_current_frame_and_time(int &step, float &t)
{
    bool bOk;
    int ret = xtc_get_current_frame_number(xtc->handle,natoms,&bOk);
    if(!bOk || ret<0) throw PterosError("Can't get current frame number");
    step = ret/xtc->steps_per_frame;
    t = xtc_get_current_frame_time(xtc->handle,natoms,&bOk);
    if(!bOk || t<0) throw PterosError("Can't get current frame time");
}

void XtcFile::tell_last_frame_and_time(int &step, float &t)
{
    step = xtc->num_frames;
    t = xtc->max_t;
}

void XtcFile::do_write(const Selection &sel, const FileContent &what)
{
    // Set box
    pteros_box_to_gmx(sel.box(),xtc->box);
    const Frame& fr = sel.get_system()->frame(sel.get_frame());
    // We need local storage, not temporary, since we pass a pointer
    auto matr = sel.get_xyz();
    rvec* x = (rvec*)matr.data(); // Pointer to storage

    int ret = write_xtc(xtc->handle,sel.size(),xtc->step,fr.time,xtc->box,x,1000);
    if(ret!=exdrOK) throw PterosError("Unable to write XTC frame {}", xtc->step);

    ++xtc->step;
}



