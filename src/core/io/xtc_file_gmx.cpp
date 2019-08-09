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


#include "xtc_file_gmx.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;
using namespace Eigen;


void XTC_file::open(char open_mode)
{
    handle = open_xtc(file_name.c_str(),&open_mode);
    if(!handle) throw Pteros_error("Unable to open TRR file {}", file_name);

    // Prepare the box just in case
    box[0][0] = 0.0; box[0][1] = 0.0; box[0][2] = 0.0;
    box[1][0] = 0.0; box[1][1] = 0.0; box[1][2] = 0.0;
    box[2][0] = 0.0; box[2][1] = 0.0; box[2][2] = 0.0;

    step = 0; // For writing
}

XTC_file::~XTC_file()
{
    close_xtc(handle);
}

bool XTC_file::do_read(System *sys, Frame *frame, const Mol_file_content &what){
    float prec;
    bool bok;

    auto ok = read_next_xtc(handle,sys->num_atoms(),&step,&frame->time, box, (rvec*)frame->coord.data() ,&prec,&bok);

    if(!bok) throw Pteros_error("XTC frame is corrupted!");

    if(ok){
        // Get box and time
        Eigen::Matrix3f b;
        for(int i=0;i<3;++i)
            for(int j=0;j<3;++j)
                b(i,j) = box[j][i];
        frame->box.set_matrix(b);
    }
    return ok;
}

void XTC_file::do_write(const Selection &sel, const Mol_file_content &what)
{
    // Set box
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            box[i][j] = sel.box().get_matrix()(j,i);

    const Frame& fr = sel.get_system()->frame(sel.get_frame());
    rvec* x = (rvec*)calloc(natoms,sizeof(*x));

    // Copy data to internal storage
    for(int i=0;i<natoms;++i){
        x[i][0] = sel.x(i);
        x[i][1] = sel.y(i);
        x[i][2] = sel.z(i);
    }

    write_xtc(handle,sel.size(),step,fr.time,box,x,1000);
    ++step;

    // Free buffer
    if(x) free(x);
}
