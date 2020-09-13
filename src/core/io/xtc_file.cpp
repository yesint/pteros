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

#include "xtc_file.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "gromacs_utils.h"
#include "xdr_utils.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


void XTC_file::open(char open_mode)
{
    bool bOk;
    handle = xdrfile_open(fname.c_str(),&open_mode);

    // Extract number of atoms
    int ok = xdr_xtc_get_natoms(handle,&natoms);
    if(!ok) throw Pteros_error("Can't read XTC number of atoms");

    // XTC file contains step number in terms of simulation steps, not saved frames
    // So we have to extract conversion factor
    int next = xtc_get_next_frame_number(handle,natoms);
    int cur = xtc_get_current_frame_number(handle,natoms,&bOk);
    if(cur<0 || next<0 || !bOk) throw Pteros_error("Can't detect number of steps per frame");
    steps_per_frame = next-cur;

    // Get total number of frames in the trajectory
    num_frames = xdr_xtc_get_last_frame_number(handle,natoms,&bOk);
    if(num_frames<0 || !bOk) throw Pteros_error("Can't get number of frames");
    num_frames /= steps_per_frame;

    // Get time step
    dt = xdr_xtc_estimate_dt(handle,natoms,&bOk);
    if(!bOk) throw Pteros_error("Can't get time step");    

    max_t = xdr_xtc_get_last_frame_time(handle,natoms,&bOk);
    if(!bOk || max_t<0) throw Pteros_error("Can't get last frame time");

    LOG()->debug("There are {} frames, max_t= {}, dt={}",num_frames,max_t,dt);




    if(!handle) throw Pteros_error("Unable to open XTC file {}", fname);

    // Prepare the box just in case
    init_gmx_box(box);

    // -1 for reading means initialization step
    step = (open_mode=='r') ? -1 : 0;
}

XTC_file::~XTC_file()
{
    if(handle) xdrfile_close(handle);
}

bool XTC_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){
    float prec;
    int ret;

    frame->coord.resize(natoms);
    ret = read_xtc(handle,natoms,&step,&frame->time,box, (rvec*)frame->coord.data(), &prec);
    if(ret == exdrENDOFFILE) return false; // End of file
    if(ret != exdrOK){
        LOG()->warn("XTC frame {} is corrupted!",step);
        return false;
    }

    if(step==0){
        LOG()->debug("Number of atoms: {}",natoms);
        LOG()->debug("XTC precision: {}",prec);
    }

    gmx_box_to_pteros(box,frame->box);
    return true;
}


void XTC_file::seek_frame(int fr)
{
    if(fr>=num_frames) throw Pteros_error("Can't seek to frame {}, there are {} frames in this file",fr,num_frames);
    int ret = xdr_xtc_seek_frame(fr*steps_per_frame,handle,natoms);
    if(ret<0) throw Pteros_error("Error seeking to frame {}",fr,num_frames);
}

void XTC_file::seek_time(float t)
{
    if(t<0 || t>max_t) throw Pteros_error("Can't seek to time {}, last time is {}",t,max_t);
    // We assume equally spaced frames in the trajectory. It's much faster
    int ret = xdr_xtc_seek_frame(ceil(t/dt)*steps_per_frame,handle,natoms);
    //int ret = xdr_xtc_seek_time(t,handle,natoms,false);
    if(ret<0) throw Pteros_error("Can't seek to time {}",t,num_frames);
}

void XTC_file::tell_current_frame_and_time(int &step, float &t)
{
    bool bOk;
    int ret = xtc_get_current_frame_number(handle,natoms,&bOk);
    if(!bOk || ret<0) throw Pteros_error("Can't get current frame number");
    step = ret/steps_per_frame;
    t = xtc_get_current_frame_time(handle,natoms,&bOk);
    if(!bOk || t<0) throw Pteros_error("Can't get current frame time");
}

void XTC_file::tell_last_frame_and_time(int &step, float &t)
{
    step = num_frames;
    t = max_t;
}

void XTC_file::do_write(const Selection &sel, const Mol_file_content &what)
{
    // Set box
    pteros_box_to_gmx(sel.box(),box);
    const Frame& fr = sel.get_system()->frame(sel.get_frame());
    // We need local storage, not temporary, since we pass a pointer
    auto matr = sel.get_xyz();
    rvec* x = (rvec*)matr.data(); // Pointer to storage

    int ret = write_xtc(handle,sel.size(),step,fr.time,box,x,1000);
    if(ret!=exdrOK) throw Pteros_error("Unable to write XTC frame {}", step);

    ++step;
}

