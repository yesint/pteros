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


#ifndef FRAME_INFO_H
#define FRAME_INFO_H

namespace pteros {

/// Information about current frame, which is passed to Consumer for analysis
/// along with the frame itself
struct Frame_info {
    /// Counts *all* frames in trajectory starting from 0
    int absolute_frame;
    /// Current time stamp
    int absolute_time;
    /// Counts only valid frames (the frames, which fall into specified time and frame range)
    /// and are sent to processing. Starts from 0.
    int valid_frame;
    /// Time of the first processed valid frame
    int first_time;
    /// Time of the last processed valid frame
    float last_time;
    /// First processed valid frame (this is an absolute value!)
    int first_frame;
    /// Last processed valid frame (this is an absolute value!)
    int last_frame;
};

}

#endif
