#ifndef __SQUARED_EXPONENTIAL
#define __SQUARED_EXPONENTIAL

#include <cuda_runtime_api.h>

__device__ void squared_exponential(int a, int b, double tau, double ell);

#endif