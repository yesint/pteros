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

#ifndef CUDA_PTEROS_ADAPTOR_H_INCLUDED
#define CUDA_PTEROS_ADAPTOR_H_INCLUDED


class GPU_Frame {
public:
    GPU_Frame(float* data, int n);
    virtual ~GPU_Frame();
    void get();
    void translate(float* shift);
private:
    // Pointer to host data
    float* host_ptr;
    // Pointer to device data
    void* dev_ptr;
    int N;
};

#endif
