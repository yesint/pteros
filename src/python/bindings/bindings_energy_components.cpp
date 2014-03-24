/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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

#include "bindings_energy_components.h"
#include "pteros/python/bindings_util.h"
#include "pteros/core/system.h"

using namespace boost::python;
using namespace pteros;

void make_bindings_Energy_components(){

    class_<Energy_components>("Energy_components", init<>())
        .def_readwrite("total",&Energy_components::total)
        .def_readwrite("lj_14",&Energy_components::lj_14)
        .def_readwrite("q_14",&Energy_components::q_14)
        .def_readwrite("lj_sr",&Energy_components::lj_sr)
        .def_readwrite("q_sr",&Energy_components::q_sr)
    ;
}
