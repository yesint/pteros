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


#include "gro_file.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/utilities.h"
#include "system_builder.h"
#include "fmt/ostream.h"
#include <cstdlib>

using namespace std;
using namespace pteros;
using namespace Eigen;

void GroFile::do_open()
{
    if(mode=='r'){
        file_handle.open(fname,std::ios::in);
        if(!file_handle) throw PterosError("Can't open GRO file '{}' for reading. Reason: {}",fname,strerror(errno));
    } else {
        file_handle.open(fname,std::ios::out);
        if(!file_handle) throw PterosError("Can't open GRO file '{}' for writing. Reason: {}",fname,strerror(errno));
    }
}

void GroFile::do_close(){
    file_handle.close();
}

bool GroFile::do_read(System *sys, Frame *frame, const FileContent &what){
    string line;
    // Skip title line
    getline(file_handle,line);

    // Read number of atoms
    getline(file_handle,line);
    natoms = str_to_int(line);

    if(what.coord()) frame->coord.resize(natoms);

    SystemBuilder builder(sys);
    Atom tmp_atom;

    // Read atoms and coordinates
    for(size_t i=0;i<natoms;++i){
        getline(file_handle,line);
        if(what.atoms()){
            tmp_atom.resid = atoi(line.substr(0,5).c_str());
            tmp_atom.resname = line.substr(5,5);
            tmp_atom.name = line.substr(10,5);
            str_trim_in_place(tmp_atom.resname);
            str_trim_in_place(tmp_atom.name);
            // Skip index field (15,5)

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
        }

        if(what.coord()){
            frame->coord[i](0) = atof(line.substr(20,8).c_str());
            frame->coord[i](1) = atof(line.substr(28,8).c_str());
            frame->coord[i](2) = atof(line.substr(36,8).c_str());
        }
    }

    //if(what.atoms()) sys->assign_resindex();

    /* Read the box
      Format: (https://manual.gromacs.org/archive/5.0.3/online/gro.html)
      v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)

      Our arrangement:
      v1(x) v2(x) v3(x)
      v1(y) v2(y) v3(y)
      v1(z) v2(z) v3(z)

      So, the sequence of reads is:
      (0,0) (1,1) (2,2) (1,0) (2,0) (0,1) (2,1) (0,2) (1,2)
    */
    if(what.coord()){
        getline(file_handle,line);

        stringstream ss(line);

        Matrix3f box;
        box.fill(0.0);
        // Read diagonal
        ss >> box(0,0),box(1,1),box(2,2);

        // Try to read non-diagonal elements. If failed we have rectangular box.
        //  If this read fails the box is not modified, so no need to take actions
        ss >> box(1,0) >> box(2,0) >> box(0,1)
           >> box(2,1) >> box(0,2) >> box(1,2);

        frame->box.set_matrix(box);
    }
    // Report success
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
        fmt::print(file_handle,
                   "{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n",
                   resid, sel.resname(i).c_str(), sel.name(i).c_str(), ind,
                   sel.x(i), sel.y(i), sel.z(i));
    }

    // Write periodic box
    if(sel.box().is_periodic()){
        auto const& b = sel.box();
        // We store box as column-vectors, while the code below hacked from VMD use row vectors,
        // so, transpose

        // Diagonal elements
        // Use same format as Gromacs for consistency, but this is free format
        fmt::print(file_handle,
                   "{:>10.4f} {:>10.4f} {:>10.4f}",
                   b.get_element(0,0), b.get_element(1,1), b.get_element(2,2));

        // Write off-diagonal only for triclinic boxes
        if(sel.box().is_triclinic()){
            // note leading space added after diagonal
            fmt::print(file_handle,
                       " {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f}",
                       b.get_element(1,0), b.get_element(2,0), b.get_element(0,1),
                       b.get_element(2,1), b.get_element(0,2), b.get_element(1,2));
        }
        // Add training newline
        fmt::print(file_handle,"\n");
    } else {
        // No box, write zero diagonal
        fmt::print(file_handle,"0.0 0.0 0.0\n");
    }
}
