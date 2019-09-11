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


#include "xtc_file.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "gromacs_utils.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


void XTC_file::open(char open_mode)
{
#ifdef USE_GROMACS
    LOG()->debug("Using gmxlib backend for XTC");
    handle = open_xtc(fname.c_str(),&open_mode);
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
