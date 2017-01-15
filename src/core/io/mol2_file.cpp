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

#include "mol2_file.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

VMDPLUGIN_EXTERN int mol2plugin_init();
VMDPLUGIN_EXTERN int mol2plugin_register(void *v, vmdplugin_register_cb cb);
VMDPLUGIN_EXTERN int mol2plugin_fini();

MOL2_file::MOL2_file(string &fname): VMD_molfile_plugin_wrapper(fname){
    plugin = molfile_plugins["mol2"];
}
