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
    cudaMallocHost((void **) &imag_h, rows*cols*sizeof(double));

    // Allocate vectors on device
    cudaMalloc((void **) &real_d,    rows*cols*sizeof(double));
    cudaMalloc((void **) &imag_d, rows*cols*sizeof(double));
}

ComplexMatrix::ComplexMatrix(int rows) {

    // Save input parameters
    this->rows       = rows;
    this->cols       = 1;

    // Allocate vectors on host
    cudaMallocHost((void **) &real_h,    rows*cols*sizeof(double));
    cudaMallocHost((void **) &imag_h, rows*cols*sizeof(double));

    // Allocate vectors on device
    cudaMalloc((void **) &real_d,    rows*cols*sizeof(double));
    cudaMalloc((void **) &imag_d, rows*cols*sizeof(double));
}

void ComplexMatrix::free() {
    
    // Free on host
    cudaFreeHost(real_h);
    cudaFreeHost(imag_h);

    // Free on device
    cudaFree(real_d);
    cudaFree(imag_d);

}

void ComplexMatrix::toHost() {

    // Send from device to host
    cudaMemcpy(real_h,    real_d,    rows * cols * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(imag_h, imag_d, rows * cols * sizeof(double), cudaMemcpyDeviceToHost);

}

void ComplexMatrix::toDevice() {

    // Send from host to device
    cudaMemcpy(real_d,    real_h,    rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(imag_d, imag_h, rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    
}

void ComplexMatrix::setHostRealValue(int r, double val) {
    real_h[r] = val;
}

void ComplexMatrix::setHostRealValue(int r, int c, double val) {
    real_h[r*cols + c] = val;
}

__device__ void ComplexMatrix::setDeviceRealValue(int r, double val) {
    real_d[r] = val;
}

__device__ void ComplexMatrix::setDeviceRealValue(int r, int c, double val) {
    real_d[r*cols + c] = val;
}

double ComplexMatrix::getHostRealValue(int r) {
    return real_h[r];
}

double ComplexMatrix::getHostRealValue(int r, int c) {
    return real_h[r*cols + c];
}

__device__ double ComplexMatrix::getDeviceRealValue(int r) {
    return real_d[r];
}

__device__ double ComplexMatrix::getDeviceRealValue(int r, int c) {
    return real_d[r*cols + c];
}

void ComplexMatrix::setHostImagValue(int r, double val) {
    imag_h[r] = val;
}

void ComplexMatrix::setHostImagValue(int r, int c, double val) {
    imag_h[r*cols + c] = val;
}

__device__ void ComplexMatrix::setDeviceImagValue(int r, double val) {
    imag_d[r] = val;
}

__device__ void ComplexMatrix::setDeviceImagValue(int r, int c, double val) {
    imag_d[r*cols + c] = val;
}

double ComplexMatrix::getHostImagValue(int r) {
    return imag_h[r];
}

double ComplexMatrix::getHostImagValue(int r, int c) {
    return imag_h[r*cols + c];
}

__device__ double ComplexMatrix::getDeviceImagValue(int r) {
    return imag_d[r];
}

__device__ double ComplexMatrix::getDeviceImagValue(int r, int c) {
    return imag_d[r*cols + c];
}

}