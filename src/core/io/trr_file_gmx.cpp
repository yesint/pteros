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


#include "trr_file_gmx.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "gromacs_utils.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


void TRR_file::open(char open_mode)
{
    handle = gmx_trr_open(fname.c_str(),&open_mode);
    if(!handle) throw Pteros_error("Unable to open TRR file {}", fname);

    // Prepare the box just in case
    init_gmx_box(box);

    step = 0; // For writing
}

TRR_file::~TRR_file()
{
    gmx_trr_close(handle);
}

bool TRR_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){

    gmx_bool bok;
    auto ok = gmx_trr_read_frame_header(handle,&header,&bok);

    if(!bok) Pteros_error("Incomplete frame header in TRR file", file_name);
    if(!ok) return false;

    rvec* x = nullptr;
    rvec* v = nullptr;
    rvec* f = nullptr;
    if(header.x_size){
        frame->coord.resize(header.natoms);
        x = (rvec*)frame->coord.data();
    }
    if(header.v_size){
        frame->vel.resize(header.natoms);
        v = (rvec*)frame->vel.data();
    }
    if(header.f_size){
        frame->force.resize(header.natoms);
        f = (rvec*)frame->force.data();
    }

    ok = gmx_trr_read_frame_data(handle,&header,box,x,v,f);

    if(ok){
        // Get box and time
        gmx_box_to_pteros(box,frame->box);
        frame->time = header.t;
    }
    return ok;
}

void TRR_file::do_write(const Selection &sel, const Mol_file_content &what)
{
    // Set box    
    pteros_box_to_gmx(sel.box(),box);

    const Frame& fr = sel.get_system()->frame(sel.get_frame());

    rvec* x = (rvec*)sel.get_xyz().data();

    rvec* v = nullptr;
    if(fr.has_vel())
        v = (rvec*)sel.get_vel().data();

    rvec* f = nullptr;
    if(fr.has_force())
        f = (rvec*)sel.get_force().data();

    gmx_trr_write_frame(handle,step,fr.time,0,box,sel.size(),x,v,f);
    ++step;
}
