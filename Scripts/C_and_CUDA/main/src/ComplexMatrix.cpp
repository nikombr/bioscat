
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../lib/ComplexMatrix.h"
using namespace std;

ComplexMatrix::ComplexMatrix(int rows, int cols) {

    // Save input parameters
    this->rows       = rows;
    this->cols       = cols;

    // Allocate vectors on host
    cudaMallocHost((void **) &real_h,    rows*cols*sizeof(double));
    cudaMallocHost((void **) &complex_h, rows*cols*sizeof(double));

    // Allocate vectors on device
    cudaMalloc((void **) &real_d,    rows*cols*sizeof(double));
    cudaMalloc((void **) &complex_d, rows*cols*sizeof(double));
}

ComplexMatrix::~ComplexMatrix() {
    
    // Free on host
    cudaFreeHost(real_h);
    cudaFreeHost(complex_h);

    // Free on device
    cudaFree(real_d);
    cudaFree(complex_d);

}

ComplexMatrix::toHost() {

    // Send from device to host
    cudaMemcpy(real_h,    real_d,    rows * cols * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(complex_h, complex_d, rows * cols * sizeof(double), cudaMemcpyDeviceToHost);

}

ComplexMatrix::toDevice() {

    // Send from host to device
    cudaMemcpy(real_d,    real_h,    rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(complex_d, complex_h, rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    
}


}