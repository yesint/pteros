#include "pteros/core/atom_proxy.h"
#include "pteros/core/system.h"
#include "pteros/core/utilities.h"

using namespace pteros;
using namespace std;

Atom_proxy::Atom_proxy(System *s, int i, int fr): atom_ptr(&s->Atom_data(i)), coord_ptr(&s->Frame_data(fr).coord[i]) {}

void Atom_proxy::set(System *s, int i, int fr){
    atom_ptr = &s->Atom_data(i);
    coord_ptr = &s->Frame_data(fr).coord[i];
}

string Atom_proxy::element_name() const {
    return get_element_name(atom_ptr->element_number);
}

float Atom_proxy::vdw() const {
    return get_vdw_radius(atom_ptr->element_number,atom_ptr->name);
}
