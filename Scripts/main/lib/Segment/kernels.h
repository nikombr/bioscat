#ifndef _KERNELS_H
#define _KERNELS_H
extern "C" {
#include "../utils/RealMatrix.h"


__global__ void setConstantKernel(RealMatrix output, int output_shift, int n, double value);

__global__ void setLinearKernel(RealMatrix output, int output_shift, int n, double value);

__global__ void setLinearKernel2(RealMatrix output, int output_shift, int n, double value);

__global__ void setReversedKernel(RealMatrix output, int output_shift, int n, double value);

__global__ void setReversedVectorKernel(RealMatrix output, int output_shift, int n, RealMatrix input, int input_shift);

__global__ void setVectorKernel(RealMatrix output, int output_shift, int n, RealMatrix input, int input_shift);

void setConstantCPU(RealMatrix output, int output_shift, int n, double value);

void setLinearCPU(RealMatrix output, int output_shift, int n, double value);

void setLinearCPU2(RealMatrix output, int output_shift, int n, double value);

void setReversedCPU(RealMatrix output, int output_shift, int n, double value);

void setReversedVectorCPU(RealMatrix output, int output_shift, int n, RealMatrix input, int input_shift);

 void setVectorCPU(RealMatrix output, int output_shift, int n, RealMatrix input, int input_shift);
}

#endif