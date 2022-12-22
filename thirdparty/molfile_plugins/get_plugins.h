#pragma once

#include "molfile_plugin.h"

extern "C" {

extern molfile_plugin_t* pdb_get_plugin_ptr();
extern molfile_plugin_t* xyz_get_plugin_ptr();
extern molfile_plugin_t* dcd_get_plugin_ptr();
#ifdef USE_TNG
extern molfile_plugin_t* tng_get_plugin_ptr();
#endif

}
