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

void Segment::computeScatteredFieldMatrices(RealMatrix x, RealMatrix y) {
    
    int rows = y.rows;
    int cols = n_int;

    ComplexMatrix F1, F2, F3;
    double const1, const2, k0;

    k0 = constants.k0;

    if (polarisation == 1) {
        F1 = E_scat_matrix.z;
        F2 = H_scat_matrix.x;
        F3 = H_scat_matrix.y;
        const1 = 1.0;
        const2 = 1/constants.eta0;
    }
    else if (polarisation == 2) {
        F1 = H_scat_matrix.z;
        F2 = E_scat_matrix.x;
        F3 = E_scat_matrix.y;
        const1 = - 1.0;
        const2 = 1/constants.eta0;
    } 
    else {
        printf("Please input 1 or 2 for the polarisation in \"computeScatteredFieldMatrices\"!\n");
        return;
    }
    if (deviceComputation) {

        // Blocks and threads
        dim3 dimBlock(32,16);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x, (cols + dimBlock.y - 1)/dimBlock.y);
        
        computeFieldMatricesKernel<<<dimGrid, dimBlock>>>(F1, F2, F3, x, aux_int.x, y, aux_int.y, const1, const2, rows, cols, k0);
        cudaCheckError();
        cudaDeviceSynchronize();
        
    }
    else {
        /*#pragma omp parallel for collapse(2) 
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                double abs_int, abs_int_ref, xdiff, ydiff, ydiff_ref, H_real, H_imag, H_real_ref, H_imag_ref, val;

                // Get data
                xdiff       = x.getHostValue(r) - x_int.getHostValue(c);
                ydiff       = y.getHostValue(r) - y_int.getHostValue(c);
                ydiff_ref   = y.getHostValue(r) + y_int.getHostValue(c);
                abs_int     = std::sqrt(xdiff * xdiff + ydiff     * ydiff);
                abs_int_ref = std::sqrt(xdiff * xdiff + ydiff_ref * ydiff_ref);

                // Compute first Hankel functions
                H_real     = H02_real(k0 * abs_int);
                H_real_ref = H02_real(k0 * abs_int_ref);
                H_imag     = H02_imag(k0 * abs_int);
                H_imag_ref = H02_imag(k0 * abs_int_ref);
                
                val = H_real + Gamma_ref * H_real_ref;
                F1.setHostRealValue(r, c, val);
                val = H_imag + Gamma_ref * H_imag_ref;
                F1.setHostImagValue(r, c, val);

                // Compute second Hankel functions
                H_real     = H12_real(k0 * abs_int);
                H_real_ref = H12_real(k0 * abs_int_ref);
                H_imag     = H12_imag(k0 * abs_int);
                H_imag_ref = H12_imag(k0 * abs_int_ref);

                val = constant * (1/abs_int      * H_imag     * ydiff + \
                        Gamma_ref * 1/abs_int_ref * H_imag_ref * ydiff_ref);
                F2.setHostRealValue(r, c, val);
                val = -constant * (1/abs_int     * H_real     * ydiff + \
                        Gamma_ref * 1/abs_int_ref * H_real_ref * ydiff_ref);
                F2.setHostImagValue(r, c, val);

                val = -constant * xdiff * (1/abs_int      * H_imag      + \
                                Gamma_ref * 1/abs_int_ref  * H_imag_ref);
                F3.setHostRealValue(r, c, val);
                val = constant * xdiff * (1/abs_int     * H_real      + \
                            Gamma_ref * 1/abs_int_ref * H_real_ref);
                F3.setHostImagValue(r, c, val);
            }
        }*/
    }
}


}