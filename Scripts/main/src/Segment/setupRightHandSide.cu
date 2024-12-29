#include <stdlib.h>
#include <stdio.h>
//#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "Segment.h"
#include "RealMatrix.h"
#include "ComplexMatrix.h"
extern "C" {
using namespace std;

void setupRightHandSide_CPU(RealMatrix b, int n_test, RealMatrix n_y, ComplexMatrix F1, ComplexMatrix F2) {

    int real_shift =  0;
    int imag_shift =  2 * n_test;
    double val;
    #pragma omp parallel for
    for (int j = 0; j < n_test; j++) {
        val =  - F1.getHostRealValue(j);
        b.setHostValue(j + real_shift, val);
        val =  - F1.getHostImagValue(j);
        b.setHostValue(j + imag_shift, val);

        val  = F2.getHostRealValue(j);
        val *= n_y.getHostValue(j);
        b.setHostValue(j + n_test + real_shift, val);
        val  = F2.getHostImagValue(j);
        val *= n_y.getHostValue(j);
        b.setHostValue(j + n_test + imag_shift, val);
    }

}

__global__ void setupRightHandSideKernel(RealMatrix b, int n_test, RealMatrix n_y, ComplexMatrix F1, ComplexMatrix F2) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n_test) {
        int real_shift =  0;
        int imag_shift =  2 * n_test;
        double val;
        val =  - F1.getDeviceRealValue(j);
        b.setDeviceValue(j + real_shift, val);
        val =  - F1.getDeviceImagValue(j);
        b.setDeviceValue(j + imag_shift, val);

        val  = F2.getDeviceRealValue(j);
        val *= n_y.getDeviceValue(j);
        b.setDeviceValue(j + n_test + real_shift, val);
        val  = F2.getDeviceImagValue(j);
        val *= n_y.getDeviceValue(j);
        b.setDeviceValue(j + n_test + imag_shift, val);
    }
}
void setupRightHandSide_GPU(RealMatrix b, int n_test, RealMatrix n_y, ComplexMatrix F1, ComplexMatrix F2) {

    // Blocks and threads
    dim3 dimBlock(256);
    dim3 dimGrid((n_test + dimBlock.x - 1)/dimBlock.x);
    
    // Run loops
    setupRightHandSideKernel<<<dimGrid, dimBlock>>>(b, n_test, n_y, F1, F2);

    // Synchronize threads before combining
    cudaDeviceSynchronize();
}


void Segment::setupRightHandSide() {
    

    if (deviceComputation) {
        if (polarisation == 1) {
            setupRightHandSide_GPU(b, n_test, normal_vectors.y, E_inc.z, H_inc.x);
        }
        else if (polarisation == 2) {
            setupRightHandSide_GPU(b, n_test, normal_vectors.y, H_inc.z, E_inc.x);
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