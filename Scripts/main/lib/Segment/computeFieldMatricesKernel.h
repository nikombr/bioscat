#ifndef _COMPUTE_FIELD_MATRICES_KERNELS_H
#define _COMPUTE_FIELD_MATRICES_KERNELS_H

extern "C" {


__host__ __device__ double H02_real(double x);
__host__ __device__ double H02_imag(double x);
__host__ __device__ double H12_real(double x);
__host__ __device__ double H12_imag(double x);

__global__ void computeFieldMatricesKernel(ComplexMatrix F1, ComplexMatrix F2, ComplexMatrix F3, RealMatrix x, RealMatrix x_aux, RealMatrix y, RealMatrix y_aux, double const1, double const2, int rows, int cols, double k);
void computeFieldMatricesCPU(ComplexMatrix F1, ComplexMatrix F2, ComplexMatrix F3, RealMatrix x, RealMatrix x_aux, RealMatrix y, RealMatrix y_aux, double const1, double const2, int rows, int cols, double k);

}

#endif