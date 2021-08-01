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

using namespace std;
using namespace pteros;
using namespace Eigen;

void GroFile::open(char open_mode)
{
    if(open_mode=='r'){
        f.open(fname.c_str(),ios_base::in);
        if(!f) throw PterosError("Can't open GRO file '{}' for reading",fname);
    } else {
        f.open(fname.c_str(),ios_base::out);
        if(!f) throw PterosError("Can't open GRO file '{}' for writing",fname);
    }
}

void GroFile::close(){
    if(f){
        f.close();
    }
}

bool GroFile::do_read(System *sys, Frame *frame, const FileContent &what){

    string line;

    int N,i,j;
    float v;
    // tmp atom
    Atom tmp_atom;
    // Temporary coordinates
    Vector3f tmp_coor;
    // Buffers for scanf
    char resname_buf[5], name_buf[5];

    // Skip header line
    getline(f,line);    

    // Read number of atoms
    getline(f,line);
    sscanf(line.c_str(),"%d",&N);

    frame->coord.resize(N);

    SystemBuilder builder(sys);

    // Read coordinates
    for(i=0;i<N;++i){
        getline(f,line);

        tmp_atom.resid = atoi(line.substr(0,5).c_str());
        tmp_atom.resname = line.substr(5,5);
        tmp_atom.name = line.substr(10,5);
        // dum - 5 chars
        tmp_coor(0) = atof(line.substr(20,8).c_str());
        tmp_coor(1) = atof(line.substr(28,8).c_str());
        tmp_coor(2) = atof(line.substr(36,8).c_str());

        str_trim_in_place(tmp_atom.resname);
        str_trim_in_place(tmp_atom.name);

        // Coordinates are in nm, so no need to convert

        if(what.atoms()){
            // Assign masses
            get_element_from_atom_name(tmp_atom.name, tmp_atom.atomic_number, tmp_atom.mass);
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
            // Add column of coordinates
            frame->coord[i] = tmp_coor;
        }
    }

    if(what.atoms()) sys->assign_resindex();

    if(what.coord()){
        // Read box. Adapted form VMD.
        stringstream ss;
        getline(f,line);
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
        // Box is in nm in gro files, no need to convert
        // Transpose the box because we want column-vectors (the code above uses row-vectors)
        box.transposeInPlace();

        frame->box.set_matrix(box);
    }

    return true;
}

void GroFile::do_write(const Selection &sel, const FileContent &what){
    int n = sel.size();
    char ch[80];

    if(!(what.coord() && what.atoms()))
        throw PterosError("It is impossible to write individual components to GRO file!");

    // Print title
    f << "Created by Pteros" << endl;
    // Write number of atoms
    f << n << endl;
    int ind,resid;
    for(int i=0;i<n;i++){
        ind = (i%99999)+1; // Prevents overflow of index field. It's not used anyway.
        resid = (sel.resid(i)%99999); // Prevents overflow of resid field.
        sprintf(ch,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f",
                //  ^Resid
                //       ^Resname
                //          ^Name
                //             ^ind
                resid, sel.resname(i).c_str(), sel.name(i).c_str(), ind,
                sel.x(i), sel.y(i), sel.z(i));
        f << ch << endl;
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
    // We are writing dimensions in nm to be compatible with Gromacs
    // Write diagonal anyway
    f << b(0,0) << " "
      << b(1,1) << " "
      << b(2,2);
    // Write off-diagonal only for triclinic boxes
    if(sel.box().is_triclinic()){
        f << " "
          << b(0,1) << " "
          << b(0,2) << " "
          << b(1,0) << " "
          << b(1,2) << " "
          << b(2,0) << " "
          << b(2,1);
    }
    // Mandatory endline at the end of file!
    f << endl;
}




