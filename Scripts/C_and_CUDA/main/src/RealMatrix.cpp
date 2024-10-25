
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../lib/RealMatrix.h"
using namespace std;

RealMatrix::RealMatrix(int rows, int cols) {

    // Save input parameters
    this->rows       = rows;
    this->cols       = cols;

    // Allocate vectors on host
    cudaMallocHost((void **) &val_h,    rows*cols*sizeof(double));

    // Allocate vectors on device
    cudaMalloc((void **) &val_d,    rows*cols*sizeof(double));

    val_h[0] = 1;
}

void RealMatrix::freeRealMatrix() {
    
    // Free on host
    cudaFreeHost(val_h);

    // Free on device
    cudaFree(val_d);

    printf("Destructed real matrices!\n");

}

void RealMatrix::toHost() {

    // Send from device to host
    cudaMemcpy(val_h,    val_d,    rows * cols * sizeof(double), cudaMemcpyDeviceToHost);

}

void RealMatrix::toDevice() {

    // Send from host to device
    cudaMemcpy(val_d,    val_h,    rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    
}


}