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

#include "pteros/python/bindings_util.h"
#include "bindings_system.h"
#include "bindings_selection.h"
#include "bindings_options_tree.h"
#include "bindings_frame_info.h"
#include "bindings_trajectory_processor.h"

/**********************
  Create python module
***********************/

// Translates Pteros_error to Python exception
void Pteros_error_translator(const pteros::Pteros_error& e) {
  PyErr_SetString(PyExc_UserWarning, const_cast<pteros::Pteros_error&>(e).what().c_str());
}


BOOST_PYTHON_MODULE(pteros)
{
    using namespace pteros;

    // Required!
    import_array();    
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

    // Register exception translator
    register_exception_translator<Pteros_error>(&Pteros_error_translator);

    // Bindings for System   
    make_bindings_System();

    // Bindings for Selection    
    make_bindings_Selection();

    // Bindings for options_tree
    make_bindings_Options_tree();

    // Bindings for Frame_info
    make_bindings_Frame_info();

    // Bindings for Trajectory_processor
    make_bindings_Trajectory_processor();

}

