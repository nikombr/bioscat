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
void Segment::computeInteriorFields() {

    int rows = n_obs;
    int cols = n_ext;
    

    if (deviceComputation) {
        //E_int_matrix.toDevice();
        //H_int_matrix.toDevice();
        double start = omp_get_wtime();
        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x);
        if (polarisation == 1) {
            computeFieldsFromMatricesKernel<<<dimGrid, dimBlock>>>(D, E_int.z, H_int.x, H_int.y, E_int_matrix.z, H_int_matrix.x, H_int_matrix.y, rows, cols);
        }
        else if (polarisation == 2) {
            computeFieldsFromMatricesKernel<<<dimGrid, dimBlock>>>(D, H_int.z, E_int.x, E_int.y, H_int_matrix.z, E_int_matrix.x, E_int_matrix.y, rows, cols);
        }
        cudaCheckError();

        

        cudaDeviceSynchronize();
        double end = omp_get_wtime();
  

    }
    else {
        if (polarisation == 1) {
            computeFieldsFromMatricesCPU(D, E_int.z, H_int.x, H_int.y, E_int_matrix.z, H_int_matrix.x, H_int_matrix.y, rows, cols);
        }
        else if (polarisation == 2) {
            computeFieldsFromMatricesCPU(D, H_int.z, E_int.x, E_int.y, H_int_matrix.z, E_int_matrix.x, E_int_matrix.y, rows, cols);
        }
        
    }
}




}