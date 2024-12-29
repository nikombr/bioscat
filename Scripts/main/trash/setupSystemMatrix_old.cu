#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "../../lib/Segment.h"
extern "C" {
#include "../../lib/utils/RealMatrix.h"
using namespace std;

#define cudaCheckError() {                                          \
    cudaError_t e = cudaGetLastError();                             \
    if (e != cudaSuccess) {                                         \
        printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
}

void setupSystemMatrix_CPU(RealMatrix A, int n_test, int n_int, int n_ext, RealMatrix n_x, RealMatrix n_y, ComplexMatrix F1_scat, ComplexMatrix F1_int,  ComplexMatrix F2_scat, ComplexMatrix F2_int, ComplexMatrix F3_scat, ComplexMatrix F3_int) {
    bool host = true; 
    bool device = false;
    RealMatrix A_real = RealMatrix(2 * n_test, n_int + n_ext, host, device);
    RealMatrix A_imag = RealMatrix(2 * n_test, n_int + n_ext, host, device);

    #pragma omp parallel for collapse(2)
    for (int r = 0; r < n_test; r++) {
        for (int c = 0; c < n_int; c++) {
            double val;
            val =  F1_scat.getHostRealValue(r, c);
            A_real.setHostValue(r, c, val);
            val =  F1_scat.getHostImagValue(r, c);
            A_imag.setHostValue(r, c, val);
        }
    }
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < n_test; r++) {
        for (int c = 0; c < n_ext; c++) {
            double val;
            val =  -F1_int.getHostRealValue(r, c);
            A_real.setHostValue(r, c + n_int, val);
            val =  -F1_int.getHostImagValue(r, c);
            A_imag.setHostValue(r, c + n_int, val);
        }
       
    }
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < n_test; r++) {
        for (int c = 0; c < n_int; c++) {
            double val;
            val  = 0;
            val += - n_y.getHostValue(r) * F2_scat.getHostRealValue(r, c);
            val +=   n_x.getHostValue(r) * F3_scat.getHostRealValue(r, c);
            A_real.setHostValue(r + n_test, c, val);
            val  = 0;
            val += - n_y.getHostValue(r) * F2_scat.getHostImagValue(r, c);
            val +=   n_x.getHostValue(r) * F3_scat.getHostImagValue(r, c);
            A_imag.setHostValue(r + n_test, c, val);
        }
    }
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < n_test; r++) {
        for (int c = 0; c < n_ext; c++) {
            double val;
            val  = 0;
            val +=   n_y.getHostValue(r) * F2_int.getHostRealValue(r, c);
            val += - n_x.getHostValue(r) * F3_int.getHostRealValue(r, c);
            A_real.setHostValue(r + n_test, c + n_int, val);
            val  = 0;
            val +=   n_y.getHostValue(r) * F2_int.getHostImagValue(r, c);
            val += - n_x.getHostValue(r) * F3_int.getHostImagValue(r, c);
            A_imag.setHostValue(r + n_test, c + n_int, val);
        }
    }

    #pragma omp parallel for collapse(2)
    for (int r = 0; r < 2 * n_test; r++) {
        for (int c = 0; c < n_ext + n_int; c++) {
            A.setHostValue(r,              c,                   A_real.getHostValue(r,c));
            A.setHostValue(r,              c + n_ext + n_int, - A_imag.getHostValue(r,c));
            A.setHostValue(r + 2 * n_test, c,                   A_imag.getHostValue(r,c));
            A.setHostValue(r + 2 * n_test, c + n_ext + n_int,   A_real.getHostValue(r,c));
        }
    }

    A_real.free();
    A_imag.free();

}

__global__ void computeScatteringAmatrix(RealMatrix A_real, RealMatrix A_imag, int n_test, int n_int, RealMatrix n_x, RealMatrix n_y, ComplexMatrix F1_scat, ComplexMatrix F2_scat, ComplexMatrix F3_scat) {
    double val;
    //printf("in kernel\n");
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    int c = threadIdx.y + blockIdx.y * blockDim.y;
    if (r < n_test && c < n_int) {
        //printf("her1\n");
        val =  F1_scat.getDeviceRealValue(r, c);
        A_real.setDeviceValue(r, c, val);
        val =  F1_scat.getDeviceImagValue(r, c);
        A_imag.setDeviceValue(r, c, val);

        val  = 0;
        val += - n_y.getDeviceValue(r) * F2_scat.getDeviceRealValue(r, c);
        val +=   n_x.getDeviceValue(r) * F3_scat.getDeviceRealValue(r, c);
        A_real.setDeviceValue(r + n_test, c, val);
        val  = 0;
        val += - n_y.getDeviceValue(r) * F2_scat.getDeviceImagValue(r, c);
        val +=   n_x.getDeviceValue(r) * F3_scat.getDeviceImagValue(r, c);
        A_imag.setDeviceValue(r + n_test, c, val);
    }
}

__global__ void computeInteriorAmatrix(RealMatrix A_real, RealMatrix A_imag, int n_test, int n_int, int n_ext, RealMatrix n_x, RealMatrix n_y, ComplexMatrix F1_int, ComplexMatrix F2_int, ComplexMatrix F3_int) {
    double val;
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    int c = threadIdx.y + blockIdx.y * blockDim.y;
    int n = n_int + n_ext;
    if (r < n_test && c < n_ext) {

        val =  -F1_int.getDeviceRealValue(r, c);
        A_real.setDeviceValue(r, c + n_int, val);
        val =  -F1_int.getDeviceImagValue(r, c);
        A_imag.setDeviceValue(r, c + n_int, val);

        val  = 0;
        val +=   n_y.getDeviceValue(r) * F2_int.getDeviceRealValue(r, c);
        val += - n_x.getDeviceValue(r) * F3_int.getDeviceRealValue(r, c);
        A_real.setDeviceValue(r + n_test, c + n_int, val);
        val  = 0;
        val +=   n_y.getDeviceValue(r) * F2_int.getDeviceImagValue(r, c);
        val += - n_x.getDeviceValue(r) * F3_int.getDeviceImagValue(r, c);
        A_imag.setDeviceValue(r + n_test, c + n_int, val);

    }
}

__global__ void combineAmatrixLoop(RealMatrix A, RealMatrix A_real, RealMatrix A_imag, int n_test, int n_int,int n_ext) {
    //printf("combine\n");
    double val;
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    int c = threadIdx.y + blockIdx.y * blockDim.y;
    int n = n_int + n_ext;
    int nA = 2*n;
    if (r < 2*n_test && c < n_int + n_ext) {
        
        A.setDeviceValue(r,              c,                   A_real.getDeviceValue(r,c));
        A.setDeviceValue(r,              c + n_ext + n_int, - A_imag.getDeviceValue(r,c));
        A.setDeviceValue(r + 2 * n_test, c,                   A_imag.getDeviceValue(r,c));
        A.setDeviceValue(r + 2 * n_test, c + n_ext + n_int,   A_real.getDeviceValue(r,c));

    }
}


void setupSystemMatrix_GPU(RealMatrix A, int n_test, int n_int, int n_ext, RealMatrix n_x, RealMatrix n_y, ComplexMatrix F1_scat, ComplexMatrix F1_int, ComplexMatrix F2_scat, ComplexMatrix F2_int, ComplexMatrix F3_scat, ComplexMatrix F3_int, RealMatrix A_real, RealMatrix A_imag) {
    bool host = false;
    bool device = true;
    //RealMatrix A_real = RealMatrix(2 * n_test, n_int + n_ext, host, device);
    //RealMatrix A_imag = RealMatrix(2 * n_test, n_int + n_ext, host, device);

    // Blocks and threads
    dim3 dimBlock(32,32);
    dim3 dimGrid((n_test + dimBlock.x - 1)/dimBlock.x, (n_int + dimBlock.y - 1)/dimBlock.y);
   
    //testker<<<dimGrid, dimBlock>>>();
    cudaDeviceSynchronize();

    computeScatteringAmatrix<<<dimGrid, dimBlock>>>(A_real, A_imag, n_test, n_int, n_x,  n_y, \
                                                    F1_scat, F2_scat, F3_scat);
    cudaCheckError(); 
    cudaDeviceSynchronize();
    // Blocks and threads
    dimGrid.y = (n_ext  + dimBlock.y - 1)/dimBlock.y;
   
    computeInteriorAmatrix<<<dimGrid, dimBlock>>>(A_real, A_imag, n_test, n_int, n_ext, n_x, n_y, \
                                     F1_int, F2_int, F3_int);
    cudaCheckError(); 
    // Synchronize threads
    cudaDeviceSynchronize();

    // Blocks and threads
    dimGrid.x = (2*n_test + dimBlock.x - 1)/dimBlock.x;
    dimGrid.y = (n_int + n_ext + dimBlock.y - 1)/dimBlock.y;
 
    combineAmatrixLoop<<< dimGrid, dimBlock>>>(A, A_real, A_imag, n_test, n_int, n_ext);
    cudaCheckError(); 
    // Synchronize threads
    cudaDeviceSynchronize();

    //A_real.free();
    //A_imag.free();

}

void Segment::setupSystemMatrix() {

    if (deviceComputation) {
        if (polarisation == 1) {
            setupSystemMatrix_GPU(A, n_test, n_int, n_ext, normal_vectors.x, normal_vectors.y, E_scat_matrix.z, E_int_matrix.z, H_scat_matrix.x, H_int_matrix.x, H_scat_matrix.y, H_int_matrix.y, A_real, A_imag);
        }
        else if (polarisation == 2) {
            setupSystemMatrix_GPU(A, n_test, n_int, n_ext, normal_vectors.x, normal_vectors.y, H_scat_matrix.z, H_int_matrix.z, E_scat_matrix.x, E_int_matrix.x, E_scat_matrix.y, E_int_matrix.y, A_real, A_imag);
        }
        else {
            printf("You have to choose either 1 or 2 as the polarisation!\n");
            return;
        }
    }
    else {
         if (polarisation == 1) {
            setupSystemMatrix_CPU(A, n_test, n_int, n_ext, normal_vectors.x, normal_vectors.y, E_scat_matrix.z, E_int_matrix.z, H_scat_matrix.x, H_int_matrix.x, H_scat_matrix.y, H_int_matrix.y);
        }
        else if (polarisation == 2) {
            setupSystemMatrix_CPU(A, n_test, n_int, n_ext, normal_vectors.x, normal_vectors.y, H_scat_matrix.z, H_int_matrix.z, E_scat_matrix.x, E_int_matrix.x, E_scat_matrix.y, E_int_matrix.y);
        }
        else {
            printf("You have to choose either 1 or 2 as the polarisation!\n");
            return;
        }

    }

    return;
}


}