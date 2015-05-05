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

#include "pteros/python/bindings_util.h"
#include "bindings_system.h"
#include "bindings_selection.h"
#include "bindings_options.h"
#include "bindings_frame_info.h"
#include "bindings_trajectory_processor.h"
#include "bindings_atom.h"
#include "bindings_atom_proxy.h"
#include "bindings_distance_search.h"
#include "bindings_frame.h"
#include "bindings_periodic_box.h"
#include "bindings_energy_components.h"
#include "bindings_utilities.h"

/**********************
  Create python module
***********************/

/*
  About bindings:
  If we write boost.python converters then one excessive copy operation
  is required in all functions, which return Eigen objects.
  Converters can not be written in other way. In principle this is not a big deal since our
  largest Eigen object is Matrix4f. However, we still need manual wrappers for some methods
  because converters can't handle non-const reference params. Also because of ambiquity
  of int and bool some methods with optional args become broken and also require
  wrappers. std::vector is also not converted, etc.
  Thus we don't use converters and wrap all Eigen and other non-trivial stuff manually.
 */

// Translates Pteros_error to Python exception
void Pteros_error_translator(const pteros::Pteros_error& e) {
  PyErr_SetString(PyExc_UserWarning, const_cast<pteros::Pteros_error&>(e).what().c_str());
}

BOOST_PYTHON_MODULE(_pteros)
{
    using namespace pteros;
    using namespace boost::python;
    using namespace boost::python::numeric;

    // Required!
    import_array();    
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    // Register exception translator
    register_exception_translator<Pteros_error>(&Pteros_error_translator);

    // Bindings for Frame
    make_bindings_Frame();

    // Bindings for Periodic_box
    make_bindings_Periodic_box();

    // Bindings for System   
    make_bindings_System();

    // Bindings for Selection    
    make_bindings_Selection();

    // Bindings for options_tree
    //make_bindings_Options_tree();

    make_bindings_Options();

    // Bindings for Frame_info
    make_bindings_Frame_info();

    // Bindings for Trajectory_processor
    make_bindings_Trajectory_processor();

    // Bindings for atom
    make_bindings_Atom();

    // Bindings for Atom_proxy
    make_bindings_Atom_proxy();

    // Bindings for grid search
    make_bindings_distance_search();

    // Energy components
    make_bindings_Energy_components();

    // Energy components
    make_bindings_utilities();
}

