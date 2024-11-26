#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
#include "../../lib/Segment/computeFieldMatricesKernel.h"
extern "C" {
using namespace std;

#define cudaCheckError() {                                          \
    cudaError_t e = cudaGetLastError();                             \
    if (e != cudaSuccess) {                                         \
        printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
}


void Segment::computeInteriorFieldMatrices(RealMatrix x, RealMatrix y) {
    
    int rows = y.rows;
    int cols = n_ext;

    ComplexMatrix F1, F2, F3;
    double const1, const2, k1;

    k1 = constants.k1;

    if (polarisation == 1) {
        F1 = E_int_matrix.z;
        F2 = H_int_matrix.x;
        F3 = H_int_matrix.y;
        const1 = 1.0;
        const2 = constants.n1/constants.eta0;
    }
    else if (polarisation == 2) {
        F1 = H_int_matrix.z;
        F2 = E_int_matrix.x;
        F3 = E_int_matrix.y;
        const1 = - 1.0;
        const2 = 1.0/(constants.eta0*constants.n1);
    } 
    else {
        printf("Please input 1 or 2 for the polarisation!\n");
    }

    if (deviceComputation) {
        // Blocks and threads
        dim3 dimBlock(32,16);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x, (cols + dimBlock.y - 1)/dimBlock.y);
        
        computeFieldMatricesKernel<<<dimGrid, dimBlock>>>(F1, F2, F3, x, aux_ext.x, y, aux_ext.y, const1, const2, rows, cols, k1);
        cudaCheckError();
        cudaDeviceSynchronize();

    }
    else {
        /*#pragma omp parallel for collapse(2) 
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                double abs_ext, xdiff, ydiff, H_real, H_imag, val;

                // Get data
                xdiff   = x.getHostValue(r) - x_ext.getHostValue(c);
                ydiff   = y.getHostValue(r) - y_ext.getHostValue(c);
                abs_ext = std::sqrt(xdiff*xdiff + ydiff*ydiff);

                // Compute first Hankel functions
                H_real = H02_real(k1 * abs_ext);
                H_imag = H02_imag(k1 * abs_ext);
                
                val = H_real;
                F1.setHostRealValue(r, c, val);
                val = H_imag;
                F1.setHostImagValue(r, c, val);

                // Compute second Hankel functions
                H_real = H12_real(k1 * abs_ext);
                H_imag = H12_imag(k1 * abs_ext);

                val =   constant * 1/abs_ext * ydiff * H_imag;
                F2.setHostRealValue(r, c, val);
                val = - constant * 1/abs_ext * ydiff * H_real;
                F2.setHostImagValue(r, c, val);

                val = -constant * 1/abs_ext * xdiff * H_imag;
                F3.setHostRealValue(r, c, val);
                val =  constant * 1/abs_ext * xdiff * H_real;
                F3.setHostImagValue(r, c, val);
            }
        }*/
    }
}


}