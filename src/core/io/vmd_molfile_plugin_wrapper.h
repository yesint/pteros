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

#ifndef VMD_MOLFILE_PLUGIN_WRAPPER_H
#define VMD_MOLFILE_PLUGIN_WRAPPER_H

#include <string>
#include "pteros/core/system.h"
#include "pteros/core/selection.h"
#include "pteros/core/mol_file.h"
#include "molfile_plugin.h"

namespace pteros {

/// Generic API for reading and writing any molecule file formats
class VMD_molfile_plugin_wrapper: public Mol_file {
public:
    // High-level API        
    VMD_molfile_plugin_wrapper(std::string& fname);
    virtual void open(char open_mode);
    virtual ~VMD_molfile_plugin_wrapper();

protected:       
    void* handle; // Handle for reading
    void* w_handle; // Handle for writing

    char mode;

    virtual bool do_read(System *sys, Frame *frame, const Mol_file_content& what);
    virtual void do_write(const Selection &sel, const Mol_file_content& what);

    molfile_plugin_t* plugin;

    int register_cb(void *v, vmdplugin_t *p);

    static std::map<std::string,molfile_plugin_t*> molfile_plugins;
};

}
#endif /* MOL_FILE_H */
