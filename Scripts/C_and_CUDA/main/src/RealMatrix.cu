
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../lib/RealMatrix.h"
using namespace std;

RealMatrix::RealMatrix(int rows) {
     // Save input parameters
    this->rows       = rows;
    this->cols       = 1;

    // Allocate vectors on host
    cudaMallocHost((void **) &val_h,    rows*cols*sizeof(double));

    // Allocate vectors on device
    cudaMalloc((void **) &val_d,    rows*cols*sizeof(double));

    val_h[0] = 1;
}

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

void RealMatrix::free() {
    
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

void RealMatrix::setHostValue(int r, double val) {
    val_h[r] = val;
}

void RealMatrix::setHostValue(int r, int c, double val) {
    val_h[r*cols + c] = val;
}

__device__ void RealMatrix::setDeviceValue(int r, double val) {
    val_h[r] = val;
}

__device__ void RealMatrix::setDeviceValue(int r, int c, double val) {
    val_h[r*cols + c] = val;
}

double RealMatrix::getHostValue(int r) {
    return val_h[r];
}

double RealMatrix::getHostValue(int r, int c) {
    return val_h[r*cols + c];
}

__device__ double RealMatrix::getDeviceValue(int r) {
    return val_h[r];
}

__device__ double RealMatrix::getDeviceValue(int r, int c) {
    return val_h[r*cols + c];
}


}