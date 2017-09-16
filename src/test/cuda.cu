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
/*
#include <stdio.h>
#include "cuda_pteros_adaptor.h"
#include <iostream>

const int N = 500000;
const int Nrun = 50000;

int main(void)
{

  float v[N*3];
  for(int i=0;i<N*3;++i) v[i]=i;

  float shift[3] = {1,2,3};

  // GPU
  GPU_Frame f(v,N);
  for(int c=0;c<Nrun;++c){
    f.translate(shift);
  }
  f.get();


 //cudaDeviceSynchronize();

  //for(int i=0;i<N*3;++i) std::cout << v[i] << " " << std::endl;

  return 0;
}
*/

#include <iostream>
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#define N 1000000
int main(void)
{
    float time_cpu;
    float time_gpu;
    int *a = new int[N];
    int *b = new int[N];
    int *c = new int[N];
    for(int i=0;i<N;i++)
    {
        a[i]=i;
        b[i]=i*i;
    }
    clock_t start_cpu,stop_cpu;
    start_cpu=clock();
    for(int i=0;i<N;i++)
    {
        c[i]=a[i]+b[i];
    }
    stop_cpu=clock();
    time_cpu=(double)(stop_cpu-start_cpu)/CLOCKS_PER_SEC;
    std::cout<<"Time to generate (CPU):"<<time_cpu<<std::endl;



    thrust::device_vector<int> X(N);
    thrust::device_vector<int> Y(N);
    thrust::device_vector<int> Z(N);
    for(int i=0;i<N;i++)
    {
        X[i]=i;
        Y[i]=i*i;
    }
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start,0);

    thrust::transform(X.begin(), X.end(),
        Y.begin(),
        Z.begin(),
        thrust::plus<int>());

    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    float elapsedTime;
    cudaEventElapsedTime(&elapsedTime,start,stop);
    std::cout<<"Time to generate (thrust):"<<elapsedTime<<std::endl;
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    return 0;
}
