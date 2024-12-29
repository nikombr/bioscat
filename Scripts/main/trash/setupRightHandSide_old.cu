#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
#include "../../lib/utils/ComplexMatrix.h"
extern "C" {
using namespace std;

void setupRightHandSide_CPU(RealMatrix b, int n_test, RealMatrix n_y, ComplexMatrix F1, ComplexMatrix F2) {
    bool host = true;
    bool device = false;
    RealMatrix b_imag = RealMatrix(2 * n_test, host, device);
    RealMatrix b_real = RealMatrix(2 * n_test, host, device);
    double val;
    #pragma omp parallel for
    for (int j = 0; j < n_test; j++) {
        val =  - F1.getHostRealValue(j);
        b_real.setHostValue(j, val);
        val =  - F1.getHostImagValue(j);
        b_imag.setHostValue(j, val);
    }
    
    #pragma omp parallel for
    for (int j = 0; j < n_test; j++) {
        val  = F2.getHostRealValue(j);
        val *= n_y.getHostValue(j);
        b_real.setHostValue(j + n_test, val);
        val  = F2.getHostImagValue(j);
        val *= n_y.getHostValue(j);
        b_imag.setHostValue(j + n_test, val);
    }

    #pragma omp parallel for
    for (int r = 0; r < 2*n_test; r++) {
        b.setHostValue(r,            b_real.getHostValue(r));
        b.setHostValue(r + 2*n_test, b_imag.getHostValue(r));
    }

    b_real.free();
    b_imag.free();
}

__global__ void firstLoop(RealMatrix b_real, int n_test, ComplexMatrix F1) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n_test) {
        double val =  - F1.getDeviceRealValue(j);
        b_real.setDeviceValue(j, val);
    }
}

__global__ void secondLoop(RealMatrix b_imag, int n_test, ComplexMatrix F1) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n_test) {
        double val =  - F1.getDeviceImagValue(j);
        b_imag.setDeviceValue(j, val);
    }
}

__global__ void thirdLoop(RealMatrix b_real, int n_test, RealMatrix n_y, ComplexMatrix F2) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n_test) {
        double val  = F2.getDeviceRealValue(j);
        val *= n_y.getDeviceValue(j);
        b_real.setDeviceValue(j + n_test, val);
    }
}

__global__ void fourthLoop(RealMatrix b_imag, int n_test, RealMatrix n_y, ComplexMatrix F2) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n_test) {
        double val  = F2.getDeviceImagValue(j);
        val *= n_y.getDeviceValue(j);
        b_imag.setDeviceValue(j + n_test, val);;
    }
}

__global__ void combineReal(RealMatrix b, RealMatrix b_real, int n_test) {
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    if (r < 2*n_test) {
        b.setDeviceValue(r, b_real.getDeviceValue(r));
    }
}

__global__ void combineImag(RealMatrix b, RealMatrix b_imag, int n_test) {
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    if (r < 2*n_test) {
        b.setDeviceValue(r + 2*n_test, b_imag.getDeviceValue(r));
    }
}

void setupRightHandSide_GPU(RealMatrix b, int n_test, RealMatrix n_y, ComplexMatrix F1, ComplexMatrix F2, RealMatrix b_real, RealMatrix b_imag) {
    //bool host = false;
    //bool device = true;
    //RealMatrix b_imag = RealMatrix(2 * n_test, host, device);
    //RealMatrix b_real = RealMatrix(2 * n_test, host, device);
    double val;

    // Blocks and threads
    dim3 dimBlock(256);
    dim3 dimGrid((n_test + dimBlock.x - 1)/dimBlock.x);
    
    // Run loops
    firstLoop<<< dimGrid, dimBlock>>>(b_real, n_test,      F1);
    secondLoop<<<dimGrid, dimBlock>>>(b_imag, n_test,      F1);
    thirdLoop<<< dimGrid, dimBlock>>>(b_real, n_test, n_y, F2);
    fourthLoop<<<dimGrid, dimBlock>>>(b_imag, n_test, n_y, F2);

    // Synchronize threads before combining
    cudaDeviceSynchronize();

    // Blocks and threads
    dimGrid.x = (2*n_test + dimBlock.x - 1)/dimBlock.x;

    // Combine results
    combineReal<<<dimGrid, dimBlock>>>(b, b_real, n_test);
    combineImag<<<dimGrid, dimBlock>>>(b, b_imag, n_test);
    // Synchronize threads
    cudaDeviceSynchronize();

    //b_real.free();
    //b_imag.free();
}


void Segment::setupRightHandSide() {
    

    if (deviceComputation) {
        if (polarisation == 1) {
            setupRightHandSide_GPU(b, n_test, normal_vectors.y, E_inc.z, H_inc.x, b_real, b_imag);
        }
        else if (polarisation == 2) {
            setupRightHandSide_GPU(b, n_test, normal_vectors.y, H_inc.z, E_inc.x, b_real, b_imag);
        }
        else {
            printf("You have to choose either 1 or 2 as the polarisation!\n");
            return;
        }
    }
    else {
        if (polarisation == 1) {
            setupRightHandSide_CPU(b, n_test, normal_vectors.y, E_inc.z, H_inc.x);
        }
        else if (polarisation == 2) {
            setupRightHandSide_CPU(b, n_test, normal_vectors.y, H_inc.z, E_inc.x);
        }
        else {
            printf("You have to choose either 1 or 2 as the polarisation!\n");
            return;
        }

        
    }
 
    
    return;
 
}

}