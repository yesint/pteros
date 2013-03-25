/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/

#include "bindings_atom.h"

void make_bindings_Atom(){

    class_<Atom>("Atom", init<>())
        .def_readwrite("resid",&Atom::resid)
        .def_readwrite("name",&Atom::name)
        .def_readwrite("chain",&Atom::chain)
        .def_readwrite("resname",&Atom::resname)
        .def_readwrite("tag",&Atom::tag)
        .def_readwrite("occupancy",&Atom::occupancy)
        .def_readwrite("beta",&Atom::beta)
        .def_readwrite("resindex",&Atom::resindex)
        .def_readwrite("mass",&Atom::mass)
        .def_readwrite("type",&Atom::type)
        .def_readwrite("type_name",&Atom::type_name)
        .def_readwrite("V",&Atom::V)
        .def_readwrite("W",&Atom::W)
    ;
}
