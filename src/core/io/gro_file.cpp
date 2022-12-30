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


#include "gro_file.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/utilities.h"
#include "system_builder.h"
#include "scn/scn.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

void GroFile::open(char open_mode)
{
    if(open_mode=='r'){
        file_handle = fopen(fname.c_str(),"r");
        if(!file_handle) throw PterosError("Can't open GRO file '{}' for reading",fname);
    } else {
        file_handle = fopen(fname.c_str(),"w");
        if(!file_handle) throw PterosError("Can't open GRO file '{}' for writing",fname);
    }
}

void GroFile::close(){
    if(file_handle){
        fclose(file_handle);
    }
}

bool GroFile::do_read(System *sys, Frame *frame, const FileContent &what){
    int N,i;
    // tmp atom
    Atom tmp_atom;
    if(what.coord()) frame->coord.resize(N);

    scn::file f(file_handle);
    auto res = scn::make_result(f);

    // Skip title line
    res = scn::ignore_until(res.range(),'\n');;
    if(!res) throw PterosError("Can't read title from GRO file! {}",
                               res.error().msg());

    // Read number of atoms
    res = scn::scan_default(res.range(),N);
    if(!res) throw PterosError("Can't read number of atoms from GRO file! {}",
                               res.error().msg());

    SystemBuilder builder(sys);

    // Read atoms and coordinates
    string dum;
    for(i=0;i<N;++i){
        if(what.atoms()){
            //fmt::print("1i={}\n",i);
            res = scn::scan(res.range(),"{:5i}{:5s}{:5s}",
                            tmp_atom.resid,
                            tmp_atom.resname,
                            tmp_atom.name);
            if(!res) throw PterosError("Corrupted atom in GRO file! Line {}: {}",
                                       i,res.error().msg());
            // Skip index field
            res = scn::ignore_until_n(res.range(),5,'\n');

            // Assign masses
            get_element_from_atom_name(tmp_atom.name,
                                       tmp_atom.atomic_number,
                                       tmp_atom.mass);
            tmp_atom.type = -1; //Undefined type so far
            // There is no chain, occupancy and beta in GRO file, so add it manually
            tmp_atom.chain = 'X';
            tmp_atom.beta = 0.0;
            tmp_atom.occupancy = 0.0;
            // We have to deduce the element number
            tmp_atom.atomic_number = get_element_number(tmp_atom.name);
            // Add new atom to the system
            builder.add_atom(tmp_atom);
        } else {
            res = scn::ignore_until_n(res.range(),20,'\n');
        }        

        if(what.coord()){
            res = scn::scan(res.range(),"{:8f}{:8f}{:8f}",
                            frame->coord[i](0),
                            frame->coord[i](1),
                            frame->coord[i](2));
            if(!res) throw PterosError("Corrupted ccordinates in GRO file! Line {}: {}",
                                       i,res.error().msg());

        }

        // Ignore until the end of line (in case of velocities)
        res = scn::ignore_until(res.range(),'\n');

        //fmt::print("2i={} {} {}\n",tmp_atom.name, tmp_atom.resname, frame->coord[i](0));
    }

    if(what.atoms()) sys->assign_resindex();

    /*
    if(what.coord()){
        // Read box. Adapted form VMD.
        stringstream ss;
        getline(file_handle,line);
        ss.clear();
        ss.str(line);
        //ss >> &x[0], &y[1], &z[2], &x[1], &x[2], &y[0], &y[2], &z[0], &z[1])
        Matrix3f box;
        box.fill(0.0);
        ss >> box(0,0) >> box(1,1) >> box(2,2);
        // Try to read nex val. If failed we have rectangular box.
        ss >> v;
        if(ss.good()){
            box(0,1) = v;
            ss >> box(0,2) >> box(1,0) >> box(1,2)
               >> box(2,0) >> box(2,1);
        }
        // Transpose the box because we want column-vectors (the code above uses row-vectors)
        box.transposeInPlace();

        frame->box.set_matrix(box);
    }
    */

    return true;
}

void GroFile::do_write(const Selection &sel, const FileContent &what){
    int n = sel.size();

    if(!(what.coord() && what.atoms()))
        throw PterosError("It is impossible to write individual components to GRO file!");

    // Print title
    fmt::print(file_handle, "Created by Pteros\n");
    // Write number of atoms
    fmt::print(file_handle,"{}\n",n);
    for(int i=0;i<n;i++){
        int ind = (i%99999)+1; // Prevents overflow of index field. It's not used anyway.
        int resid = (sel.resid(i)%99999); // Prevents overflow of resid field.
        fmt::print(file_handle,"{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n",
                   resid, sel.resname(i).c_str(), sel.name(i).c_str(), ind,
                   sel.x(i), sel.y(i), sel.z(i));
    }

    // Write periodic box
    Eigen::Matrix3f b;
    if(sel.box().is_periodic()){
        // We store box as column-vectors, while the code below hacked from VMD use row vectors,
        // so, transpose
        b = sel.box().get_matrix().transpose();
    } else {
        b.fill(0.0);
    }
    /*
    // We are writing dimensions in nm to be compatible with Gromacs
    // Write diagonal anyway
    file_handle << b(0,0) << " "
      << b(1,1) << " "
      << b(2,2);
    // Write off-diagonal only for triclinic boxes
    if(sel.box().is_triclinic()){
        file_handle << " "
          << b(0,1) << " "
          << b(0,2) << " "
          << b(1,0) << " "
          << b(1,2) << " "
          << b(2,0) << " "
          << b(2,1);
    }
    */
    // Mandatory endline at the end of file!
    //file_handle << endl;
}




