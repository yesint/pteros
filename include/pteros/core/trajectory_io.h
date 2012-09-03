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

#ifndef TRAJECTORY_IO_H
#define TRAJECTORY_IO_H

#include "pteros/core/generic_trajectory.h"
#include "pteros/core/format_recognition.h"

namespace pteros {

/// Factory returning trajectory reader for given file format
boost::shared_ptr<Generic_trajectory> trajectory_reader_factory(FILE_FORMATS fmt);

}

#endif // TRAJECTORY_IO_H
