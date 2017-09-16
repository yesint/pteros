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

#include "cuda_pteros_adaptor.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>

GPU_Frame::GPU_Frame(float *data, int n)
{
    // Save host ptr
    host_ptr = data;
    // Allocate device memory
    cudaMalloc(&dev_ptr, n*sizeof(float)*3);
    // Remember size
    N = n;
    // Copy host to device
    cudaMemcpy(dev_ptr, (void*)host_ptr, N*sizeof(float)*3, cudaMemcpyHostToDevice);
}

GPU_Frame::~GPU_Frame()
{
    if(dev_ptr) cudaFree(dev_ptr);
}

void GPU_Frame::get()
{
    // Copy back to host
    cudaMemcpy((void*)host_ptr, dev_ptr, N*sizeof(float)*3, cudaMemcpyDeviceToHost);
}

struct do_shift{
    float s[3];
    do_shift(float* shift){
        s[0] = shift[0];
        s[1] = shift[1];
        s[2] = shift[2];
    }
    __device__ float3 operator()(float3& v){
        float3 res;
        res.x = v.x+s[0];
        res.y = v.y+s[1];
        res.z = v.z+s[2];
        return res;
    }
};

void GPU_Frame::translate(float *shift)
{
    // Map thrust vector to device memory
    thrust::device_ptr<float3> p((float3*)dev_ptr);
    // Do operation
    thrust::transform(p,p+N,p,do_shift(shift));
}
