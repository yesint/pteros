#include "pteros/core/atom_handler.h"
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/utilities.h"
#include "pteros/core/pteros_error.h"

using namespace pteros;
using namespace std;

AtomHandler::AtomHandler(const Selection &sel, int i){
    set(sel,i);
}

AtomHandler::AtomHandler(const System &sys, int i, int fr){
    set(sys,i,fr);
}

void AtomHandler::set(const Selection &sel, int i)
{
    if(sel.get_frame()>=sel.get_system()->num_frames())
        throw PterosError("Can't point at frame {}, there are only {} frames in the system!",i,sel.get_system()->num_frames());
    coord_ptr = sel.get_system()->frame(sel.get_frame()).coord.begin();
    coord_ptr += i;
    atom_ptr = sel.get_system()->atoms_begin();
    atom_ptr += i;
    ind = i;
}

void AtomHandler::set(const System &sys, int i, int fr)
{
    if(fr>=sys.num_frames())
        throw PterosError("Can't point ot frame {}, there are only {} frames in the system!",i,sys.num_frames());
    coord_ptr = const_cast<System&>(sys).frame(fr).coord.begin();
    coord_ptr += i;
    atom_ptr = const_cast<System&>(sys).atoms_begin();
    atom_ptr += i;
    ind = i;
}

string AtomHandler::element_name() const {
    return get_element_name(atom_ptr->atomic_number);
}

float AtomHandler::vdw() const {
    return get_vdw_radius(atom_ptr->atomic_number,atom_ptr->name);
}
