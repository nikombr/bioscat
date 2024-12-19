#include <stdlib.h>
#include <stdio.h>
//#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
#include "../../lib/Segment/computeFieldsFromMatricesKernel.h"
extern "C" {
using namespace std;

#define cudaCheckError() {                                          \
    cudaError_t e = cudaGetLastError();                             \
    if (e != cudaSuccess) {                                         \
        printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
}
void Segment::computeScatteredFields() {

    int rows = n_obs;
    int cols = n_int;

    if (deviceComputation) {
  
        double start = omp_get_wtime();
        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x);
        if (polarisation == 1) {
            computeFieldsFromMatricesKernel<<<dimGrid, dimBlock>>>(C, E_scat.z, H_scat.x, H_scat.y, E_scat_matrix.z, H_scat_matrix.x, H_scat_matrix.y, rows, cols);
        }
        else if (polarisation == 2) {
            computeFieldsFromMatricesKernel<<<dimGrid, dimBlock>>>(C, H_scat.z, E_scat.x, E_scat.y, H_scat_matrix.z, E_scat_matrix.x, E_scat_matrix.y, rows, cols);
        }
        cudaCheckError();
        cudaDeviceSynchronize();
        double end = omp_get_wtime();

    }
    else {
        if (polarisation == 1) {
            computeFieldsFromMatricesCPU(C, E_scat.z, H_scat.x, H_scat.y, E_scat_matrix.z, H_scat_matrix.x, H_scat_matrix.y, rows, cols);
        }
        else if (polarisation == 2) {
            computeFieldsFromMatricesCPU(C, H_scat.z, E_scat.x, E_scat.y, H_scat_matrix.z, E_scat_matrix.x, E_scat_matrix.y, rows, cols);
        }
        
    }
}




}