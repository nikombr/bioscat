#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "Segment.h"
extern "C" {
#include "RealMatrix.h"
using namespace std;

#define cudaCheckError() {                                          \
    cudaError_t e = cudaGetLastError();                             \
    if (e != cudaSuccess) {                                         \
        printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
}

void setupSystemMatrix_CPU(RealMatrix A, int n_test, int n_int, int n_ext, RealMatrix n_x, RealMatrix n_y, ComplexMatrix F1_scat, ComplexMatrix F1_int,  ComplexMatrix F2_scat, ComplexMatrix F2_int, ComplexMatrix F3_scat, ComplexMatrix F3_int) {

    int real_rshift_1 = 0;
    int real_cshift_1 = 0;

    int real_rshift_2 = 2 * n_test;
    int real_cshift_2 = n_ext + n_int;

    int imag_rshift_1 = 0;
    int imag_cshift_1 = n_ext + n_int;

    int imag_rshift_2 = 2 * n_test;
    int imag_cshift_2 = 0;


    #pragma omp parallel for collapse(2)
    for (int r = 0; r < n_test; r++) {
        for (int c = 0; c < n_int; c++) {
            double val;
            val =  F1_scat.getHostRealValue(r, c);
            A.setHostValue(r + real_rshift_1, c + real_cshift_1, val);
            A.setHostValue(r + real_rshift_2, c + real_cshift_2, val);
            val =  F1_scat.getHostImagValue(r, c);
            A.setHostValue(r + imag_rshift_1, c + imag_cshift_1, - val);
            A.setHostValue(r + imag_rshift_2, c + imag_cshift_2, val);

            val  = 0;
            val += - n_y.getHostValue(r) * F2_scat.getHostRealValue(r, c);
            val +=   n_x.getHostValue(r) * F3_scat.getHostRealValue(r, c);
            A.setHostValue(r + n_test + real_rshift_1, c + real_cshift_1, val);
            A.setHostValue(r + n_test + real_rshift_2, c + real_cshift_2, val);
            val  = 0;
            val += - n_y.getHostValue(r) * F2_scat.getHostImagValue(r, c);
            val +=   n_x.getHostValue(r) * F3_scat.getHostImagValue(r, c);
            A.setHostValue(r + n_test + imag_rshift_1, c + imag_cshift_1, - val);
            A.setHostValue(r + n_test + imag_rshift_2, c + imag_cshift_2, val);
        }
    }

    #pragma omp parallel for collapse(2)
    for (int r = 0; r < n_test; r++) {
        for (int c = 0; c < n_ext; c++) {
            double val;
            val =  -F1_int.getHostRealValue(r, c);
            A.setHostValue(r + real_rshift_1, c + n_int + real_cshift_1, val);
            A.setHostValue(r + real_rshift_2, c + n_int + real_cshift_2, val);
            val =  -F1_int.getHostImagValue(r, c);
            A.setHostValue(r + imag_rshift_1, c + n_int + imag_cshift_1, - val);
            A.setHostValue(r + imag_rshift_2, c + n_int + imag_cshift_2, val);

            val  = 0;
            val +=   n_y.getHostValue(r) * F2_int.getHostRealValue(r, c);
            val += - n_x.getHostValue(r) * F3_int.getHostRealValue(r, c);
            A.setHostValue(r + n_test + real_rshift_1, c + n_int + real_cshift_1, val);
            A.setHostValue(r + n_test + real_rshift_2, c + n_int + real_cshift_2, val);
            val  = 0;
            val +=   n_y.getHostValue(r) * F2_int.getHostImagValue(r, c);
            val += - n_x.getHostValue(r) * F3_int.getHostImagValue(r, c);
            A.setHostValue(r + n_test + imag_rshift_1, c + n_int + imag_cshift_1, - val);
            A.setHostValue(r + n_test + imag_rshift_2, c + n_int + imag_cshift_2, val);
        }
       
    }
   
}

__global__ void computeScatteringAmatrix(RealMatrix A, int n_test, int n_int, int n_ext, RealMatrix n_x, RealMatrix n_y, ComplexMatrix F1_scat, ComplexMatrix F2_scat, ComplexMatrix F3_scat) {
    double val;
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    int c = threadIdx.y + blockIdx.y * blockDim.y;
    if (r < n_test && c < n_int) {
        int real_rshift_1 = 0;
        int real_cshift_1 = 0;

        int real_rshift_2 = 2 * n_test;
        int real_cshift_2 = n_ext + n_int;

        int imag_rshift_1 = 0;
        int imag_cshift_1 = n_ext + n_int;

        int imag_rshift_2 = 2 * n_test;
        int imag_cshift_2 = 0;
        
        val =  F1_scat.getDeviceRealValue(r, c);
        A.setDeviceValue(r + real_rshift_1, c + real_cshift_1, val);
        A.setDeviceValue(r + real_rshift_2, c + real_cshift_2, val);
        val =  F1_scat.getDeviceImagValue(r, c);
        A.setDeviceValue(r + imag_rshift_1, c + imag_cshift_1, - val);
        A.setDeviceValue(r + imag_rshift_2, c + imag_cshift_2, val);

        val  = 0;
        val += - n_y.getDeviceValue(r) * F2_scat.getDeviceRealValue(r, c);
        val +=   n_x.getDeviceValue(r) * F3_scat.getDeviceRealValue(r, c);
        A.setDeviceValue(r + n_test + real_rshift_1, c + real_cshift_1, val);
        A.setDeviceValue(r + n_test + real_rshift_2, c + real_cshift_2, val);
        val  = 0;
        val += - n_y.getDeviceValue(r) * F2_scat.getDeviceImagValue(r, c);
        val +=   n_x.getDeviceValue(r) * F3_scat.getDeviceImagValue(r, c);
        A.setDeviceValue(r + n_test + imag_rshift_1, c + imag_cshift_1, - val);
        A.setDeviceValue(r + n_test + imag_rshift_2, c + imag_cshift_2, val);
    }
}

__global__ void computeInteriorAmatrix(RealMatrix A, int n_test, int n_int, int n_ext, RealMatrix n_x, RealMatrix n_y, ComplexMatrix F1_int, ComplexMatrix F2_int, ComplexMatrix F3_int) {
    double val;
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    int c = threadIdx.y + blockIdx.y * blockDim.y;
    int n = n_int + n_ext;
    if (r < n_test && c < n_ext) {
        int real_rshift_1 = 0;
        int real_cshift_1 = 0;

        int real_rshift_2 = 2 * n_test;
        int real_cshift_2 = n_ext + n_int;

        int imag_rshift_1 = 0;
        int imag_cshift_1 = n_ext + n_int;

        int imag_rshift_2 = 2 * n_test;
        int imag_cshift_2 = 0;

        val =  -F1_int.getDeviceRealValue(r, c);
        A.setDeviceValue(r + real_rshift_1, c + n_int + real_cshift_1, val);
        A.setDeviceValue(r + real_rshift_2, c + n_int + real_cshift_2, val);
        val =  -F1_int.getDeviceImagValue(r, c);
        A.setDeviceValue(r + imag_rshift_1, c + n_int + imag_cshift_1, - val);
        A.setDeviceValue(r + imag_rshift_2, c + n_int + imag_cshift_2, val);

        val  = 0;
        val +=   n_y.getDeviceValue(r) * F2_int.getDeviceRealValue(r, c);
        val += - n_x.getDeviceValue(r) * F3_int.getDeviceRealValue(r, c);
        A.setDeviceValue(r + n_test + real_rshift_1, c + n_int + real_cshift_1, val);
        A.setDeviceValue(r + n_test + real_rshift_2, c + n_int + real_cshift_2, val);
        val  = 0;
        val +=   n_y.getDeviceValue(r) * F2_int.getDeviceImagValue(r, c);
        val += - n_x.getDeviceValue(r) * F3_int.getDeviceImagValue(r, c);
        A.setDeviceValue(r + n_test + imag_rshift_1, c + n_int + imag_cshift_1, - val);
        A.setDeviceValue(r + n_test + imag_rshift_2, c + n_int + imag_cshift_2, val);

    }
}


void setupSystemMatrix_GPU(RealMatrix A, int n_test, int n_int, int n_ext, RealMatrix n_x, RealMatrix n_y, ComplexMatrix F1_scat, ComplexMatrix F1_int, ComplexMatrix F2_scat, ComplexMatrix F2_int, ComplexMatrix F3_scat, ComplexMatrix F3_int) {

    // Blocks and threads
    dim3 dimBlock(32,32);
    dim3 dimGrid((n_test + dimBlock.x - 1)/dimBlock.x, (n_int + dimBlock.y - 1)/dimBlock.y);

    computeScatteringAmatrix<<<dimGrid, dimBlock>>>(A, n_test, n_int, n_ext, n_x,  n_y, \
                                                    F1_scat, F2_scat, F3_scat);
    cudaCheckError(); 
   
    computeInteriorAmatrix<<<dimGrid, dimBlock>>>(A, n_test, n_int, n_ext, n_x, n_y, \
                                     F1_int, F2_int, F3_int);
    cudaCheckError(); 
    // Synchronize threads
    cudaDeviceSynchronize();
}

void Segment::setupSystemMatrix() {

    if (deviceComputation) {
        if (polarisation == 1) {
            setupSystemMatrix_GPU(A, n_test, n_int, n_ext, normal_vectors.x, normal_vectors.y, E_scat_matrix.z, E_int_matrix.z, H_scat_matrix.x, H_int_matrix.x, H_scat_matrix.y, H_int_matrix.y);
        }
        else if (polarisation == 2) {
            setupSystemMatrix_GPU(A, n_test, n_int, n_ext, normal_vectors.x, normal_vectors.y, H_scat_matrix.z, H_int_matrix.z, E_scat_matrix.x, E_int_matrix.x, E_scat_matrix.y, E_int_matrix.y);
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