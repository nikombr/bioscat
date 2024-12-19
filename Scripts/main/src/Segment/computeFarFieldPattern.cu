#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
extern "C" {
using namespace std;

#define cudaCheckError() {                                          \
    cudaError_t e = cudaGetLastError();                             \
    if (e != cudaSuccess) {                                         \
        printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
}

__global__ void computeFarFieldPatternKernel(ComplexMatrix F, ComplexMatrix C, RealMatrix phi, RealMatrix x_int, RealMatrix y_int, int n_obs, int n_int, double k) {
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    if (r < n_obs) {
        double exp_real, exp_imag, rho_int, phi_int, x, y, C_real, C_imag, phiIdx, val_real,val_imag,scaling_constant;
        val_real = 0;
        val_imag = 0;
        phiIdx = phi.getDeviceValue(r);
        for (int c = 0; c < n_int; c++) {

            // Get data
            x = x_int.getDeviceValue(c);
            y = y_int.getDeviceValue(c);
            
            phi_int = atan2(y, x);
            rho_int = sqrt(x * x + y * y);
            
            C_real = C.getDeviceRealValue(c);
            C_imag = C.getDeviceImagValue(c);

            // Compute first Hankel functions
            exp_real = cos(k * rho_int * cos(phiIdx - phi_int));
            exp_imag = sin(k * rho_int * cos(phiIdx - phi_int));

            // Compute the first field
            val_real += C_real * exp_real - C_imag * exp_imag;
            val_imag += C_real * exp_imag + C_imag * exp_real;
        }
        // Setup scaling constant
        scaling_constant = sqrt(2/(M_PI*k));
        F.setDeviceRealValue(r, scaling_constant*val_real);
        F.setDeviceImagValue(r, scaling_constant*val_imag);
    }
}

void Segment::computeFarFieldPattern(RealMatrix phi) {

    //int rows = phi.rows;
    //int cols = y_int.rows;
    double k0 = constants.k0;

    if (deviceComputation) {

        // Blocks and threads
        dim3 dimBlock(32);
        dim3 dimGrid((n_obs + dimBlock.x - 1)/dimBlock.x);
        computeFarFieldPatternKernel<<<dimGrid, dimBlock>>>(F, C, phi, aux_int.x, aux_int.y, n_obs, n_int, k0);
        cudaCheckError();
        cudaDeviceSynchronize();
        
    }
    else {
        for (int r = 0; r < n_obs; r++) {
            double exp_real, exp_imag, rho_int, phi_int, x, y, C_real, C_imag, phiIdx, val_real,val_imag,scaling_constant;
            val_real = 0;
            val_imag = 0;
            phiIdx = phi.getHostValue(r);
            for (int c = 0; c < n_int; c++) {

                // Get data
                x = aux_int.x.getHostValue(c);
                y = aux_int.y.getHostValue(c);
                
                phi_int = atan2(y, x);
                rho_int = sqrt(x * x + y * y);
                
                C_real = C.getHostRealValue(c);
                C_imag = C.getHostImagValue(c);

                // Compute first Hankel functions
                exp_real = cos(k0 * rho_int * cos(phiIdx - phi_int));
                exp_imag = sin(k0 * rho_int * cos(phiIdx - phi_int));

                // Compute the first field
                val_real += C_real * exp_real - C_imag * exp_imag;
                val_imag += C_real * exp_imag + C_imag * exp_real;
            }
            // Setup scaling constant
            scaling_constant = sqrt(2/(M_PI*k0));
            F.setHostRealValue(r, scaling_constant*val_real);
            F.setHostImagValue(r, scaling_constant*val_imag);
        }
    }

}

}