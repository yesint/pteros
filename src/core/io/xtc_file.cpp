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

#ifdef USE_GROMACS
#include "gromacs/fileio/xdrf.h"
#include "gromacs/fileio/xdr_datatype.h"
#include "gromacs/fileio/gmxfio_impl.h"
#include "gromacs/utility/futil.h"
#endif

using namespace std;
using namespace pteros;
using namespace Eigen;

#ifdef USE_GROMACS
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifndef XTC_MAGIC
#    define XTC_MAGIC 1995
#endif
#define XDR_INT_SIZE 4
const int header_size = 16;

static int xtc_at_header_start(FILE* fp, XDR* xdrs, int natoms, int* timestep, float* time)
{
    int       i_inp[3];
    float     f_inp[10];
    int       i;
    gmx_off_t off;


    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }
    /* read magic natoms and timestep */
    for (i = 0; i < 3; i++)
    {
        if (!xdr_int(xdrs, &(i_inp[i])))
        {
            gmx_fseek(fp, off + XDR_INT_SIZE, SEEK_SET);
            return -1;
        }
    }
    /* quick return */
    if (i_inp[0] != XTC_MAGIC)
    {
        if (gmx_fseek(fp, off + XDR_INT_SIZE, SEEK_SET))
        {
            return -1;
        }
        return 0;
    }
    /* read time and box */
    for (i = 0; i < 10; i++)
    {
        if (!xdr_float(xdrs, &(f_inp[i])))
        {
            gmx_fseek(fp, off + XDR_INT_SIZE, SEEK_SET);
            return -1;
        }
    }
    /* Make a rigourous check to see if we are in the beggining of a header
       Hopefully there are no ambiguous cases */
    /* This check makes use of the fact that the box matrix has 3 zeroes on the upper
       right triangle and that the first element must be nonzero unless the entire matrix is zero
     */
    if (i_inp[1] == natoms
        && ((f_inp[1] != 0 && f_inp[6] == 0) || (f_inp[1] == 0 && f_inp[5] == 0 && f_inp[9] == 0)))
    {
        if (gmx_fseek(fp, off + XDR_INT_SIZE, SEEK_SET))
        {
            return -1;
        }
        *time     = f_inp[0];
        *timestep = i_inp[2];
        return 1;
    }
    if (gmx_fseek(fp, off + XDR_INT_SIZE, SEEK_SET))
    {
        return -1;
    }
    return 0;
}

static int xtc_get_current_frame_number(FILE* fp, XDR* xdrs, int natoms, gmx_bool* bOK)
{
    gmx_off_t off;
    int       ret;
    int       step;
    float     time;
    *bOK = false;

    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }


    while (true)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = true;
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                *bOK = false;
                return -1;
            }
            return step;
        }
        else if (ret == -1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return -1;
        }
        else if (ret == 0)
        {
            /*Go back.*/
            if (gmx_fseek(fp, -2 * XDR_INT_SIZE, SEEK_CUR))
            {
                return -1;
            }
        }
    }
}


int xtc_get_next_frame_number(FILE* fp, XDR* xdrs, int natoms)
{
    gmx_off_t off;
    int       step;
    float     time;
    int       ret;

    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }

    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdrs, &step);
    while (true)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return step;
        }
        else if (ret == -1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
        }
    }
}

static float xtc_get_current_frame_time(FILE* fp, XDR* xdrs, int natoms, gmx_bool* bOK)
{
    gmx_off_t off;
    int       step;
    float     time;
    int       ret;
    *bOK = false;

    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }

    while (true)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = true;
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                *bOK = false;
                return -1;
            }
            return time;
        }
        else if (ret == -1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return -1;
        }
        else if (ret == 0)
        {
            /*Go back.*/
            if (gmx_fseek(fp, -2 * XDR_INT_SIZE, SEEK_CUR))
            {
                return -1;
            }
        }
    }
}

static float xtc_get_next_frame_time(FILE* fp, XDR* xdrs, int natoms, gmx_bool* bOK)
{
    gmx_off_t off;
    float     time;
    int       step;
    int       ret;
    *bOK = false;

    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }
    /* read one int just to make sure we dont read this frame but the next */
    xdr_int(xdrs, &step);
    while (true)
    {
        ret = xtc_at_header_start(fp, xdrs, natoms, &step, &time);
        if (ret == 1)
        {
            *bOK = true;
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                *bOK = false;
                return -1;
            }
            return time;
        }
        else if (ret == -1)
        {
            if (gmx_fseek(fp, off, SEEK_SET))
            {
                return -1;
            }
            return -1;
        }
    }
}

static float xdr_xtc_estimate_dt(FILE* fp, XDR* xdrs, int natoms, gmx_bool* bOK)
{
    float     res;
    float     tinit;
    gmx_off_t off;

    *bOK = false;
    if ((off = gmx_ftell(fp)) < 0)
    {
        return -1;
    }

    tinit = xtc_get_current_frame_time(fp, xdrs, natoms, bOK);

    if (!(*bOK))
    {
        return -1;
    }

    res = xtc_get_next_frame_time(fp, xdrs, natoms, bOK);

    if (!(*bOK))
    {
        return -1;
    }

    res -= tinit;
    if (0 != gmx_fseek(fp, off, SEEK_SET))
    {
        *bOK = false;
        return -1;
    }
    return res;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#endif


void XTC_file::open(char open_mode)
{
#ifdef USE_GROMACS
    LOG()->debug("Using gmxlib backend for XTC");
    handle = open_xtc(fname.c_str(),&open_mode);

    // Extract number of atoms
    int magic;
    gmx_bool bOk;
    gmx_off_t pos = gmx_ftell(handle->fp); // save pos
    bOk = xdr_int(handle->xdr,&magic);
    if(!bOk) throw Pteros_error("Can't read XTC magic number");
    bOk = xdr_int(handle->xdr,&natoms);
    if(!bOk) throw Pteros_error("Can't read XTC number of atoms");
    gmx_fseek(handle->fp,pos,SEEK_SET); // Rewind

    // XTC file contains step number in terms of simulation steps, not saved frames
    // So we have to extract conversion factor
    int next = xtc_get_next_frame_number(handle->fp,handle->xdr,natoms);
    int cur = xtc_get_current_frame_number(handle->fp,handle->xdr,natoms,&bOk);
    if(cur<0 || next<0 || !bOk) throw Pteros_error("Can't detect number of steps per frame");
    steps_per_frame = next-cur;

    // Get total number of frames in the trajectory
    num_frames = xdr_xtc_get_last_frame_number(handle->fp,handle->xdr,natoms,&bOk);
    if(num_frames<0 || !bOk) throw Pteros_error("Can't get number of frames");
    num_frames /= steps_per_frame;

    // Get time step
    dt = xdr_xtc_estimate_dt(handle->fp,handle->xdr,natoms,&bOk);
    if(!bOk) throw Pteros_error("Can't get time step");    

    max_t = xdr_xtc_get_last_frame_time(handle->fp,handle->xdr,natoms,&bOk);
    if(!bOk || max_t<0) throw Pteros_error("Can't get last frame time");

    LOG()->debug("There are {} frames, max_t= {}, dt={}",num_frames,max_t,dt);

#else
    LOG()->debug("Using xdrfile backend for XTC");
    handle = xdrfile_open(fname.c_str(),&open_mode);
#endif

    if(!handle) throw Pteros_error("Unable to open XTC file {}", fname);

    // Prepare the box just in case
    init_gmx_box(box);

    // -1 for reading means initialization step
    step = (open_mode=='r') ? -1 : 0;
}

XTC_file::~XTC_file()
{
#ifdef USE_GROMACS
    if(handle) close_xtc(handle);
#else
    if(handle) xdrfile_close(handle);
#endif

}

bool XTC_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){
    float prec;

#ifdef USE_GROMACS
    gmx_bool bok,ok;
    if(step<0){
        // First read allocates storage, so we are obliged to copy afterwards
        rvec* x;
        ok = read_first_xtc(handle,&natoms,&step,&frame->time, box, &x ,&prec,&bok);
        if(!ok) return false; // End of file
        frame->coord.resize(natoms);
        for(int i=0;i<natoms;++i)
            frame->coord[i] = Eigen::Map<Eigen::Vector3f>(x[i]);        
    } else {
        // Next reads can use frame storage directly
        frame->coord.resize(natoms);
        ok = read_next_xtc(handle,natoms,&step,&frame->time, box, (rvec*)frame->coord.data() ,&prec,&bok);
        if(!ok) return false; // End of file
    }
    if(!bok){
        LOG()->warn("XTC frame {} is corrupted!",step);
        return false;
    }
#else
    int ret;
    if(step<0){
        ret = read_xtc_natoms(const_cast<char*>(fname.c_str()), &natoms);
        if(ret!=exdrOK) throw Pteros_error("Can't read natoms from XTC");
    }

    frame->coord.resize(natoms);
    ret = read_xtc(handle,natoms,&step,&frame->time,box, (rvec*)frame->coord.data(), &prec);
    if(ret == exdrENDOFFILE) return false; // End of file
    if(ret != exdrOK){
        LOG()->warn("XTC frame {} is corrupted!",step);
        return false;
    }
#endif

    if(step==0){
        LOG()->debug("Number of atoms: {}",natoms);
        LOG()->debug("XTC precision: {}",prec);
    }

    gmx_box_to_pteros(box,frame->box);
    return true;
}


#ifdef USE_GROMACS
void XTC_file::seek_frame(int fr)
{
    if(fr>=num_frames) throw Pteros_error("Can't seek to frame {}, there are {} frames in this file",fr,num_frames);
    int ret = xdr_xtc_seek_frame(fr*steps_per_frame,handle->fp,handle->xdr,natoms);
    if(ret<0) throw Pteros_error("Error seeking to frame {}",fr,num_frames);
}

void XTC_file::seek_time(float t)
{
    if(t<0 || t>max_t) throw Pteros_error("Can't seek to time {}, last time is {}",t,max_t);
    // For some reason seek_time crashes if one didn't call current_time first. Weird...
    gmx_bool bOk;
    float tc = xtc_get_current_frame_time(handle->fp,handle->xdr,natoms,&bOk);
    if(!bOk || t<0) throw Pteros_error("Can't get current frame time");

    int ret = xdr_xtc_seek_time(t,handle->fp,handle->xdr,natoms,bOk);
    if(!bOk || ret<0) throw Pteros_error("Can't seek to time {}",t,num_frames);
}

void XTC_file::tell_current_frame_and_time(int &step, float &t)
{
    gmx_bool bOk;
    int ret = xtc_get_current_frame_number(handle->fp,handle->xdr,natoms,&bOk);
    if(!bOk || ret<0) throw Pteros_error("Can't get current frame number");
    step = ret/steps_per_frame;
    t = xtc_get_current_frame_time(handle->fp,handle->xdr,natoms,&bOk);
    if(!bOk || t<0) throw Pteros_error("Can't get current frame time");
}

void XTC_file::tell_last_frame_and_time(int &step, float &t)
{
    step = num_frames;
    t = max_t;
}
#endif

void XTC_file::do_write(const Selection &sel, const Mol_file_content &what)
{
    // Set box
    pteros_box_to_gmx(sel.box(),box);
    const Frame& fr = sel.get_system()->frame(sel.get_frame());
    // We need local storage, not temporary, since we pass a pointer
    auto matr = sel.get_xyz();
    rvec* x = (rvec*)matr.data(); // Pointer to storage

    // Name and signature of this function is the same in gmxlib and xdrfile but return codes differs
#ifdef USE_GROMACS
    int ret = write_xtc(handle,sel.size(),step,fr.time,box,x,1000);
    if(!ret) throw Pteros_error("Unable to write XTC frame {}",step);
#else
    int ret = write_xtc(handle,sel.size(),step,fr.time,box,x,1000);
    if(ret!=exdrOK) throw Pteros_error("Unable to write XTC frame {}", step);
#endif

    ++step;
}

