#include <stdlib.h>
#include <cublas_v2.h>
#include <stdio.h>
#include <math.h>
#include <curand_kernel.h>

extern "C" {
#include "../../lib/GaussianProcess.h"
#include <cblas.h>
#include <ctime>
#define M_PI 3.14159265358979323846


void GaussianProcess::generate_random_vector() {

    double U1, U2;

    for (int i = 0; i < n; i+=2) {
        U1 = ((double) rand())/((double) RAND_MAX);
        U2 = ((double) rand())/((double) RAND_MAX);
        p_h[i]   = sqrt(-2 * log(U1)) * cos(2 * M_PI * U2);
        p_h[i+1] = sqrt(-2 * log(U1)) * sin(2 * M_PI * U2);

    }

    if (device) {
        // Send to device
        cudaMemcpy(p_d, p_h, n * sizeof(double), cudaMemcpyHostToDevice);
    }


}

__global__ void scaleDown(double * p_d, int n) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < n) {
        p_d[i] = p_d[i]*1e-8;
    }
}

void GaussianProcess::realisation() {

    generate_random_vector();

    if (device) {

        status = cublasDtrmv(handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, M_log, n, p_d, 1);
        
        // Check if cublasDtrmv was successful
        if (status != CUBLAS_STATUS_SUCCESS) {
            printf("cublasDtrmv failed with error code: %d\n", status);
        }

        //dim3 dimBlock(32);
        //dim3 dimGrid((n + dimBlock.x - 1)/dimBlock.x);
        //scaleDown<<<dimGrid, dimBlock>>>(p_d, n);

        // Send to host
        //cudaMemcpy(p_h, p_d, n * sizeof(double), cudaMemcpyDeviceToHost);

        cudaDeviceSynchronize();
    }
    else {
        cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, *M_h, n, p_h, 1);
    }

}


}