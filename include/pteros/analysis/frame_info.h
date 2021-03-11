/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#pragma once

namespace pteros {

/// Information about current frame, which is passed to Consumer for analysis
/// along with the frame itself
struct FrameInfo {
    /// Counts *all* frames in trajectory starting from 0
    int absolute_frame;
    /// Current time stamp
    float absolute_time;
    /// Counts only valid frames (the frames, which fall into specified time and frame range)
    /// and are sent to processing. Starts from 0.
    int valid_frame;
    /// Time of the first processed valid frame
    float first_time;
    /// Time of the last processed valid frame
    float last_time;
    /// First processed valid frame (this is an absolute value!)
    int first_frame;
    /// Last processed valid frame (this is an absolute value!)
    int last_frame;
};

}


