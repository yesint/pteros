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


#include "trr_file.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "xdr_utils.h"
#include "xdrfile_trr.h"
#include "gmx_box_utils.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

// Internal TRR data
struct TrrFile::TrrData {
    TrrData(): handle(nullptr) {}
    ~TrrData(){
        if(handle) xdrfile_close(handle);
    }

    XDRFILE* handle;
    matrix box;
    int step;
};

TrrFile::TrrFile(const string &fname, char open_mode):
    FileHandler(fname,open_mode),
    trr(new TrrData)
{}

TrrFile::~TrrFile()
{}

void TrrFile::do_open()
{    
    trr->handle = xdrfile_open(fname.c_str(),&mode);
    if(!trr->handle) throw PterosError("Unable to open TRR file {}", fname);
    // -1 for reading means initialization step
    trr->step = (mode=='r') ? -1 : 0;
}

void TrrFile::do_close()
{}

bool TrrFile::do_read(System *sys, Frame *frame, const FileContent &what){
    bool has_x=true, has_v=false, has_f=false, ok;

    if(trr->step<0){
        // Read header only once on first step
        int xsz,vsz,fsz;
        check_trr_content(trr->handle,&natoms,&xsz,&vsz,&fsz);
        has_x = (xsz>0);
        has_v = (vsz>0);
        has_f = (fsz>0);
        if(!has_x) throw PterosError("Pteros can't read TRR files without coordinates!");
        LOG()->debug("TRR file has: x({}), v({}), f({})",has_x,has_v,has_f);
    }

    rvec* x = nullptr;
    rvec* v = nullptr;
    rvec* f = nullptr;
    if(has_x){
        frame->coord.resize(natoms);
        x = (rvec*)frame->coord.data();
    } else {
        throw PterosError("Pteros can't read TRR files without coordinates!");
    }

    if(has_v){
        frame->vel.resize(natoms);
        v = (rvec*)frame->vel.data();
    }
    if(has_f){
        frame->force.resize(natoms);
        f = (rvec*)frame->force.data();
    }

    float lambda;
    ok = read_trr(trr->handle,natoms,&trr->step,&frame->time,&lambda,trr->box,x,v,f);

    // Get box
    if(ok) gmx_box_to_pteros(trr->box,frame->box);
    return ok;
}

void TrrFile::do_write(const Selection &sel, const FileContent &what)
{
    // Set box    
    pteros_box_to_gmx(sel.box(),trr->box);

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

    int ret = write_trr(trr->handle,sel.size(),trr->step,fr.time,0,trr->box,x,v,f);
    if(ret!=exdrOK) throw PterosError("Unable to write TRR frame");

    ++trr->step;
}
