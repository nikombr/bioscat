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

void Segment::computeScatteredFieldMatrices(RealMatrix x, RealMatrix y) {

    //E_scat_matrix.setDeviceZero();
    //H_scat_matrix.setDeviceZero();
    
    int rows = y.rows;
    int cols = n_int;

    ComplexMatrix F1, F2, F3;
    double const1, const2, k0;

    k0 = constants.k0;

    if (polarisation == 1) {
        F1 = E_scat_matrix.z;
        F2 = H_scat_matrix.x;
        F3 = H_scat_matrix.y;
        const1 = 1.0;
        const2 = 1/constants.eta0;
    }
    else if (polarisation == 2) {
        F1 = H_scat_matrix.z;
        F2 = E_scat_matrix.x;
        F3 = E_scat_matrix.y;
        const1 = 1/constants.eta0;
        const2 = - 1.0;
    } 
    else {
        printf("Please input 1 or 2 for the polarisation in \"computeScatteredFieldMatrices\"!\n");
        return;
    }
    if (deviceComputation) {

        // Blocks and threads
        dim3 dimBlock(32,16);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x, (cols + dimBlock.y - 1)/dimBlock.y);
        
        computeFieldMatricesKernel<<<dimGrid, dimBlock>>>(F1, F2, F3, x, aux_int.x, y, aux_int.y, const1, const2, rows, cols, k0);
        cudaCheckError();
        cudaDeviceSynchronize();
        
    }
    else {
        computeFieldMatricesCPU(F1, F2, F3, x, aux_int.x, y, aux_int.y, const1, const2, rows, cols, k0);
    }
}


}