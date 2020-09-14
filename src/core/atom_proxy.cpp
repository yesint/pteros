#include "pteros/core/atom_proxy.h"
#include "pteros/core/system.h"
#include "pteros/core/utilities.h"

using namespace pteros;
using namespace std;

Atom_proxy::Atom_proxy(System *s, int i, int fr): atom_ptr(&s->atom(i)), coord_ptr(&s->frame(fr).coord[i]), ind(i) {
    v_ptr = (s->frame(fr).has_vel()) ? &s->frame(fr).vel[i] : nullptr;
    f_ptr = (s->frame(fr).has_force()) ? &s->frame(fr).force[i] : nullptr;
}

void Atom_proxy::set(System *s, int i, int fr){
    ind = i;
    atom_ptr = &s->atom(i);
    coord_ptr = &s->frame(fr).coord[i];
    v_ptr = (s->frame(fr).has_vel()) ? &s->frame(fr).vel[i] : nullptr;
    f_ptr = (s->frame(fr).has_force()) ? &s->frame(fr).force[i] : nullptr;
}

string Atom_proxy::element_name() const {
    return get_element_name(atom_ptr->atomic_number);
}

float Atom_proxy::vdw() const {
    return get_vdw_radius(atom_ptr->atomic_number,atom_ptr->name);
}
