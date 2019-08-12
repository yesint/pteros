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


#include "trr_file.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "gromacs_utils.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


void TRR_file::open(char open_mode)
{
#ifdef USE_GROMACS
    LOG()->debug("Using gmxlib backend for TRR");
    handle = gmx_trr_open(fname.c_str(),&open_mode);
#else
    LOG()->debug("Using xdrfile backend for TRR");
    handle = xdrfile_open(fname.c_str(),&open_mode);
#endif

    if(!handle) throw Pteros_error("Unable to open TRR file {}", fname);

    // Prepare the box just in case
    init_gmx_box(box);

    // -1 for reading means initialization step
    step = (open_mode=='r') ? -1 : 0;
}

TRR_file::~TRR_file()
{
#ifdef USE_GROMACS
    if(handle) gmx_trr_close(handle);
#else
    if(handle) xdrfile_close(handle);
#endif
}

#ifndef USE_GROMACS
// For xdrfile backend we add function to xdr source to check v and f, so extern it here
extern "C" {
    extern int check_trr_content(char *fn, int* natoms, int* xsz, int* vsz, int* fsz);
}
#endif

bool TRR_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){
    bool has_x, has_v, has_f, ok;

#ifdef USE_GROMACS
    // For gmxlib we read header for each frame
    gmx_bool bok;
    gmx_trr_header_t header;
    ok = gmx_trr_read_frame_header(handle,&header,&bok);
    if(!bok) Pteros_error("Incomplete frame header in TRR file", fname);
    if(!ok) return false;
    if(step<0){
        natoms = header.natoms;
        has_x = (header.x_size>0);
        has_v = (header.v_size>0);
        has_f = (header.f_size>0);        
    }
#else
    if(step<0){
        // For gmxlib we read header directly only once
        int xsz,vsz,fsz;
        check_trr_content(const_cast<char*>(fname.c_str()),&natoms,&xsz,&vsz,&fsz);
        has_x = (xsz>0);
        has_v = (vsz>0);
        has_f = (fsz>0);
        if(!has_x) throw Pteros_error("Pteros can't read TRR files without coordinates!");
    }
#endif

    if(step<0) LOG()->debug("TRR file has: x({}), v({}), f({})",has_x,has_v,has_f);

    rvec* x = nullptr;
    rvec* v = nullptr;
    rvec* f = nullptr;
    if(has_x){
        frame->coord.resize(natoms);
        x = (rvec*)frame->coord.data();
    } else {
        throw Pteros_error("Pteros can't read TRR files without coordinates!");
    }

    if(has_v){
        frame->vel.resize(natoms);
        v = (rvec*)frame->vel.data();
    }
    if(has_f){
        frame->force.resize(natoms);
        f = (rvec*)frame->force.data();
    }

#ifdef USE_GROMACS
    ok = gmx_trr_read_frame_data(handle,&header,box,x,v,f);
    if(ok){
        frame->time = header.t;
        step = header.step;
    }
#else
    float lambda;
    ok = read_trr(handle,natoms,&step,&frame->time,&lambda,box,x,v,f);
#endif
    // Get box
    if(ok) gmx_box_to_pteros(box,frame->box);
    return ok;
}

void TRR_file::do_write(const Selection &sel, const Mol_file_content &what)
{
    // Set box    
    pteros_box_to_gmx(sel.box(),box);

    const Frame& fr = sel.get_system()->frame(sel.get_frame());

    // We need local storage, not temporaries, since we pass a pointer
    Eigen::MatrixXf matr_x, matr_v, matr_f;

    matr_x = sel.get_xyz();
    rvec* x = (rvec*)matr_x.data(); // Pointer to storage

    rvec* v = nullptr;
    if(fr.has_vel()){
        matr_v = sel.get_vel();
        v = (rvec*)matr_v.data();
    }

    rvec* f = nullptr;
    if(fr.has_force()){
        matr_f = sel.get_force();
        f = (rvec*)matr_f.data();
    }

#ifdef USE_GROMACS
    gmx_trr_write_frame(handle,step,fr.time,0,box,sel.size(),x,v,f);
#else
    int ret = write_trr(handle,sel.size(),step,fr.time,0,box,x,v,f);
    if(ret!=exdrOK) throw Pteros_error("Unable to write TRR frame");
#endif

    ++step;
}
