#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
#include "../../lib/Segment/computeFieldMatricesKernel.h"
extern "C" {
using namespace std;

#define cudaCheckError() {                                          \
    cudaError_t e = cudaGetLastError();                             \
    if (e != cudaSuccess) {                                         \
        printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
}


void Segment::computeInteriorFieldMatrices(RealMatrix x, RealMatrix y) {
    
    int rows = y.rows;
    int cols = n_ext;

    ComplexMatrix F1, F2, F3;
    double const1, const2, k1;

    k1 = constants.k1;

    if (polarisation == 1) {
        F1 = E_int_matrix.z;
        F2 = H_int_matrix.x;
        F3 = H_int_matrix.y;
        const1 = 1.0;
        const2 = constants.n1/constants.eta0;
    }
    else if (polarisation == 2) {
        F1 = H_int_matrix.z;
        F2 = E_int_matrix.x;
        F3 = E_int_matrix.y;
        const1 = constants.n1/constants.eta0;
        const2 = - 1.0;
    } 
    else {
        printf("Please input 1 or 2 for the polarisation!\n");
    }

    if (deviceComputation) {
        // Blocks and threads
        dim3 dimBlock(32,16);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x, (cols + dimBlock.y - 1)/dimBlock.y);
        
        computeFieldMatricesKernel<<<dimGrid, dimBlock>>>(F1, F2, F3, x, aux_ext.x, y, aux_ext.y, const1, const2, rows, cols, k1);
        cudaCheckError();
        cudaDeviceSynchronize();

    }
    else {
        //#pragma omp parallel for collapse(2) 
        computeFieldMatricesCPU(F1, F2, F3, x, aux_ext.x, y, aux_ext.y, const1, const2, rows, cols, k1);
    }
}


}