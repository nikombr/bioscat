#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include "Segment.h"
#include "RealMatrix.h"
extern "C" {
using namespace std;


__global__ void computeFieldsFromMatricesKernel(ComplexMatrix C, ComplexMatrix F1, ComplexMatrix F2, ComplexMatrix F3, ComplexMatrix F1_matrix, ComplexMatrix F2_matrix, ComplexMatrix F3_matrix, int rows, int cols) {
    int k = threadIdx.x + blockIdx.x * blockDim.x;
    if (k < rows) {
        
        double f1_real, f1_imag, f2_real, f2_imag, f3_real, f3_imag;
        f1_real = 0.0; f2_real = 0.0; f3_real = 0.0;
        f1_imag = 0.0; f2_imag = 0.0; f3_imag = 0.0;
        
        for (int j = 0; j < cols; j++) {
            // Get factor
            double C_real = C.getDeviceRealValue(j);
            double C_imag = C.getDeviceImagValue(j);

            // Computing real values
            f1_real += F1_matrix.getDeviceRealValue(k,j) * C_real;
            f1_real -= F1_matrix.getDeviceImagValue(k,j) * C_imag;
            f2_real += F2_matrix.getDeviceRealValue(k,j) * C_real;
            f2_real -= F2_matrix.getDeviceImagValue(k,j) * C_imag;
            f3_real += F3_matrix.getDeviceRealValue(k,j) * C_real;
            f3_real -= F3_matrix.getDeviceImagValue(k,j) * C_imag;

            // Computing imagninary values
            f1_imag += F1_matrix.getDeviceRealValue(k,j) * C_imag;
            f1_imag += F1_matrix.getDeviceImagValue(k,j) * C_real;
            f2_imag += F2_matrix.getDeviceRealValue(k,j) * C_imag;
            f2_imag += F2_matrix.getDeviceImagValue(k,j) * C_real;
            f3_imag += F3_matrix.getDeviceRealValue(k,j) * C_imag;
            f3_imag += F3_matrix.getDeviceImagValue(k,j) * C_real;
    
        }

        F1.setDeviceRealValue(k, f1_real);
        F2.setDeviceRealValue(k, f2_real);
        F3.setDeviceRealValue(k, f3_real);
        F1.setDeviceImagValue(k, f1_imag);
        F2.setDeviceImagValue(k, f2_imag);
        F3.setDeviceImagValue(k, f3_imag);

        //if (k == 0) printf("E_scat = %f\n",C.getDeviceRealValue(0));
    }
    
}

void computeFieldsFromMatricesCPU(ComplexMatrix C, ComplexMatrix F1, ComplexMatrix F2, ComplexMatrix F3, ComplexMatrix F1_matrix, ComplexMatrix F2_matrix, ComplexMatrix F3_matrix, int rows, int cols) {
    #pragma omp parallel for
    for (int k = 0; k < rows; k++) {
        
        double f1_real, f1_imag, f2_real, f2_imag, f3_real, f3_imag;
        f1_real = 0.0; f2_real = 0.0; f3_real = 0.0;
        f1_imag = 0.0; f2_imag = 0.0; f3_imag = 0.0;
        
        for (int j = 0; j < cols; j++) {
            // Get factor
            double C_real = C.getHostRealValue(j);
            double C_imag = C.getHostImagValue(j);

            // Computing real values
            f1_real += F1_matrix.getHostRealValue(k,j) * C_real;
            f1_real -= F1_matrix.getHostImagValue(k,j) * C_imag;
            f2_real += F2_matrix.getHostRealValue(k,j) * C_real;
            f2_real -= F2_matrix.getHostImagValue(k,j) * C_imag;
            f3_real += F3_matrix.getHostRealValue(k,j) * C_real;
            f3_real -= F3_matrix.getHostImagValue(k,j) * C_imag;

            // Computing imagninary values
            f1_imag += F1_matrix.getHostRealValue(k,j) * C_imag;
            f1_imag += F1_matrix.getHostImagValue(k,j) * C_real;
            f2_imag += F2_matrix.getHostRealValue(k,j) * C_imag;
            f2_imag += F2_matrix.getHostImagValue(k,j) * C_real;
            f3_imag += F3_matrix.getHostRealValue(k,j) * C_imag;
            f3_imag += F3_matrix.getHostImagValue(k,j) * C_real;
    
        }

        F1.setHostRealValue(k, f1_real);
        F2.setHostRealValue(k, f2_real);
        F3.setHostRealValue(k, f3_real);
        F1.setHostImagValue(k, f1_imag);
        F2.setHostImagValue(k, f2_imag);
        F3.setHostImagValue(k, f3_imag);

        //if (k == 0) printf("E_scat = %f\n",C.getDeviceRealValue(0));
    }
    
}


}