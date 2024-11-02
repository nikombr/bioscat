#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
#include "../../lib/ComplexMatrix.h"
using namespace std;

void setupRightHandSide_CPU(RealMatrix b, int n_test, RealMatrix n_y, ComplexMatrix *firstField_inc, ComplexMatrix* firstField_ref, ComplexMatrix*secondField_inc, ComplexMatrix*secondField_ref) {

    RealMatrix b_imag = RealMatrix(2 * n_test);
    RealMatrix b_real = RealMatrix(2 * n_test);
    double val;
    for (int j = 0; j < n_test; j++) {
        val =  - firstField_inc->getHostRealValue(j) - firstField_ref->getHostRealValue(j);
        b_real.setHostValue(j, val);
        val =  - firstField_inc->getHostImagValue(j) - firstField_ref->getHostImagValue(j);
        b_imag.setHostValue(j, val);
    }
    
    for (int j = 0; j < n_test; j++) {
        val  = secondField_inc->getHostRealValue(j) + secondField_ref->getHostRealValue(j);
        val *= n_y.getHostValue(j);
        b_real.setHostValue(j + n_test, val);
        val  = secondField_inc->getHostImagValue(j) + secondField_ref->getHostImagValue(j);
        val *= n_y.getHostValue(j);
        b_imag.setHostValue(j + n_test, val);
    }

    for (int r = 0; r < 2*n_test; r++) {
        b.setHostValue(r,            b_real.getHostValue(r));
        b.setHostValue(r + 2*n_test, b_imag.getHostValue(r));
    }

    b_real.free();
    b_imag.free();
}

__global__ void firstLoop(RealMatrix b_real, int n_test, ComplexMatrix *firstField_inc, ComplexMatrix *firstField_ref) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n_test) {
        double val =  - firstField_inc->getDeviceRealValue(j) - firstField_ref->getDeviceRealValue(j);
        b_real.setDeviceValue(j, val);
    }
}

__global__ void secondLoop(RealMatrix b_imag, int n_test, ComplexMatrix *firstField_inc, ComplexMatrix *firstField_ref) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n_test) {
        double val =  - firstField_inc->getDeviceImagValue(j) - firstField_ref->getDeviceImagValue(j);
        b_imag.setDeviceValue(j, val);
    }
}

__global__ void thirdLoop(RealMatrix b_real, int n_test, RealMatrix n_y, ComplexMatrix *secondField_inc, ComplexMatrix *secondField_ref) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n_test) {
        double val  = secondField_inc->getDeviceRealValue(j) + secondField_ref->getDeviceRealValue(j);
        val *= n_y.getDeviceValue(j);
        b_real.setDeviceValue(j + n_test, val);
    }
}

__global__ void fourthLoop(RealMatrix b_imag, int n_test, RealMatrix n_y, ComplexMatrix *secondField_inc, ComplexMatrix *secondField_ref) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n_test) {
        double val  = secondField_inc->getDeviceImagValue(j) + secondField_ref->getDeviceImagValue(j);
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

void setupRightHandSide_GPU(RealMatrix b, int n_test, RealMatrix n_y, ComplexMatrix *firstField_inc, ComplexMatrix* firstField_ref, ComplexMatrix*secondField_inc, ComplexMatrix*secondField_ref) {
    RealMatrix b_imag = RealMatrix(2 * n_test);
    RealMatrix b_real = RealMatrix(2 * n_test);
    double val;

    // Blocks and threads
    dim3 dimBlock(256);
    dim3 dimGrid((n_test + dimBlock.x - 1)/dimBlock.x);

    // Run loops
    firstLoop<<< dimGrid, dimBlock>>>(b_real, n_test,      firstField_inc,  firstField_ref);
    secondLoop<<<dimGrid, dimBlock>>>(b_imag, n_test,      firstField_inc,  firstField_ref);
    thirdLoop<<< dimGrid, dimBlock>>>(b_real, n_test, n_y, secondField_inc, secondField_ref);
    fourthLoop<<<dimGrid, dimBlock>>>(b_imag, n_test, n_y, secondField_inc, secondField_ref);

    // Synchronize threads before combining
    cudaDeviceSynchronize();

    // Blocks and threads
    dimGrid.x = (2*n_test + dimBlock.x - 1)/dimBlock.x;

    // Combine results
    combineReal<<<dimGrid, dimBlock>>>(b, b_real, n_test);
    combineImag<<<dimGrid, dimBlock>>>(b, b_imag, n_test);

    // Synchronize threads
    cudaDeviceSynchronize();

    b_real.free();
    b_imag.free();
}


void Segment::setupRightHandSide() {
    
    
    
    // Placeholder for fields
    ComplexMatrix *firstField_inc, * firstField_ref, *secondField_inc,*secondField_ref;

    // Choose non-zero fields
    if (polarisation == 1) {
        firstField_inc = &E_inc_vector.z;
        firstField_ref = &E_ref_vector.z;
        secondField_inc = &H_inc_vector.x;
        secondField_ref = &H_ref_vector.x;
    }
    else if (polarisation == 2) {
        firstField_inc = &H_inc_vector.x;
        firstField_ref = &H_ref_vector.x;
        secondField_inc = &E_inc_vector.z;
        secondField_ref = &E_ref_vector.z;
    }
    else {
        printf("You have to choose either 1 or 2 as the polarisation!\n");
        return;
    }

    if (deviceComputation) {
        setupRightHandSide_GPU(b, n_test, n_y, firstField_inc, firstField_ref, secondField_inc, secondField_ref);
    }
    else {
        setupRightHandSide_CPU(b, n_test, n_y, firstField_inc, firstField_ref, secondField_inc, secondField_ref);
    }

    (*firstField_inc).free();
    (*firstField_ref).free();
    (*secondField_inc).free();
    (*secondField_ref).free();
    
    
    /*shift = n_test;
    for (int j = 0; j < n_top; j++) {
        val  = secondField_inc->getHostRealValue(j) + secondField_ref->getHostRealValue(j);
        val *= n_y.getHostValue(j);
        b_real.setHostValue(j + shift, val);
        val  = secondField_inc->getHostImagValue(j) + secondField_ref->getHostImagValue(j);
        val *= n_y.getHostValue(j);
        b_imag.setHostValue(j + shift, val);
    }

    shift += n_top;
    for(int j = 0; j < n_right; j++) {
        val = 0.0;
        b_real.setHostValue(j + shift, val);
        b_imag.setHostValue(j + shift, val);
    }

    shift += n_right;
    for(int j = 0; j < n_bottom; j++) {
        val  = secondField_inc->getHostRealValue(j + n_top + n_right) + secondField_ref->getHostRealValue(j + n_top + n_right);
        b_real.setHostValue(j + shift, val);
        val  = secondField_inc->getHostImagValue(j + n_top + n_right) + secondField_ref->getHostImagValue(j + n_top + n_right);
        b_imag.setHostValue(j + shift, val);
    }

    shift += n_bottom;
    for(int j = 0; j < n_left; j++) {
        val = 0.0;
        b_real.setHostValue(j + shift, val);
        b_imag.setHostValue(j + shift, val);
    }*/
    /*printf("b:\n");
    for (int j = 0; j < 2*n_test; j++) {
        printf("%e\t + i(%e)\n",b_real.getHostValue(j),b_imag.getHostValue(j));
    }*/

    
    if (n_test < 200) {
        char * filename;
        FILE *file;
    /*filename = "b_real_C.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < 2*n_test; r++) {
        fprintf(file, "%e\n", b_real.getHostValue(r));
    }
    fclose(file);
    filename = "b_imag_C.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < 2*n_test; r++) {
        fprintf(file, "%e\n", b_imag.getHostValue(r));
    }
    fclose(file);*/
    filename = "bbig_C.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < b.rows; r++) {
        fprintf(file, "%e\n", b.getHostValue(r));
    }
    fclose(file);
    }
    
    
    return;
 
}

}