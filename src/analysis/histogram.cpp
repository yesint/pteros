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

#include "pteros/analysis/histogram.h"
#include "pteros/core/pteros_error.h"

using namespace std;
using namespace pteros;


Histogram::Histogram(HISTOGRAM_TYPE t, float min_val, float max_val, int N):
    hist_type(t), maxval(max_val), minval(min_val), Nbins(N)
{
    data.resize(Nbins);
    delta = (max_val-min_val)/float(Nbins);
}
