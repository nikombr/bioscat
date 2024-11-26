#ifndef _COMPUTE_FIELDS_FROM_MATRICES_KERNELS_H
#define _COMPUTE_FIELDS_FROM_MATRICES_KERNELS_H
extern "C" {

__global__ void computeFieldsFromMatricesKernel(ComplexMatrix C, ComplexMatrix F1, ComplexMatrix F2, ComplexMatrix F3, ComplexMatrix F1_matrix, ComplexMatrix F2_matrix, ComplexMatrix F3_matrix, int rows, int cols);

}

#endif