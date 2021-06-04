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


#include "pteros/core/file_handler.h"
#include "pteros/core/pteros_error.h"

#include "pdb_file.h"
#include "dcd_file.h"
#include "gro_file.h"
#include "xyz_file.h"
#include "trr_file.h"
#include "xtc_file.h"

#ifdef USE_TNGIO
#include "tng_file.h"
#endif

#ifdef USE_OPENBABEL
#include "pdbqt_file.h"
#include "mol2_file.h"
#endif

#ifdef USE_GROMACS
#include "tpr_file.h"
#endif

using namespace std;
using namespace pteros;

FileHandler::FileHandler(string& file_name){
    fname = file_name;
}

FileHandler::~FileHandler(){
    close();
}

bool FileHandler::read(System *sys, Frame *frame, const FileContent &what){
    sanity_check_read(sys,frame,what);    
    return do_read(sys,frame,what);
}

void FileHandler::write(const Selection &sel, const FileContent &what) {
    sanity_check_write(sel,what);
    do_write(sel,what);
}

void FileHandler::sanity_check_read(System *sys, Frame *frame, const FileContent& what) const {
    auto c = get_content_type();
    if( !c.atoms() && what.atoms() )
        throw PterosError("Can't read structure from this file type!");
    if( !c.coord() && what.coord() )
        throw PterosError("Can't read coordinates from this file type!");
    if( !c.traj() && what.traj() )
        throw PterosError("Can't read coordinates from this file type!");
    if( !c.top() && what.top() )
        throw PterosError("Can't read topology from this file type!");
    if( !what.atoms() && !what.coord() && !what.traj() && !what.top() )
        throw PterosError("Nothing to read!");
    if(what.atoms() && !sys)
        throw PterosError("System should be provided to read structure!");
    if(what.top() && !sys)
        throw PterosError("System should be provided to read topology!");
    if((what.coord() || what.traj()) && !frame)
        throw PterosError("Frame should be provided to read coordinates!");
    if(what.atoms() && sys && sys->num_atoms()>0)
        throw PterosError("Can't read structure because system already has atoms!");
}

void FileHandler::sanity_check_write(const Selection &sel, const FileContent &what) const{
    if(!what.atoms() && !what.coord() && !what.traj() && !what.top())
        throw PterosError("Nothing to write!");
    if( !get_content_type().atoms() && what.atoms() )
        throw PterosError("Can't write structure from this file type!");
    if( !get_content_type().coord() && what.coord() )
        throw PterosError("Can't write coordinates from this file type!");
    if( !get_content_type().traj() && what.traj() )
        throw PterosError("Can't write coordinates from this file type!");
    if( !get_content_type().top() && what.top() )
        throw PterosError("Can't write topology from this file type!");
}

unique_ptr<FileHandler> FileHandler::recognize(string fname){
    std::string ext = fname.substr(fname.find_last_of(".") + 1);

         if(ext=="xtc")     return FileHandler_ptr(new XtcFile(fname));
    else if(ext=="trr")     return FileHandler_ptr(new TrrFile(fname));
    else if(ext=="pdb")     return FileHandler_ptr(new PdbFile(fname));
    else if(ext=="gro")     return FileHandler_ptr(new GroFile(fname));
    else if(ext=="dcd")     return FileHandler_ptr(new DcdFile(fname));
    else if(ext=="xyz")     return FileHandler_ptr(new XyzFile(fname));
#ifdef USE_TNGIO
    else if(ext=="tng")     return FileHandler_ptr(new TngFile(fname));
#endif
#ifdef USE_GROMACS
    else if(ext=="tpr")     return FileHandler_ptr(new TprFile(fname));
#endif
#ifdef USE_OPENBABEL
    else if(ext=="mol2")    return FileHandler_ptr(new Mol2File(fname));
    else if(ext=="pdbqt")   return FileHandler_ptr(new PdbqtFile(fname));
#endif
    else throw PterosError("File extension '{}' is not recognized!",ext);
}

FileHandler_ptr FileHandler::open(string fname, char open_mode)
{
    auto handle = recognize(fname);
    handle->open(open_mode);
    return handle;
}

void FileHandler::close()
{

}
