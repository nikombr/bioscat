#include <stdlib.h>
#include <stdio.h>
//#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
#include "../../lib/Segment/computeFieldsFromMatricesKernel.h"
extern "C" {
using namespace std;

#define cudaCheckError() {                                          \
    cudaError_t e = cudaGetLastError();                             \
    if (e != cudaSuccess) {                                         \
        printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
}
void Segment::computeInteriorFields() {

    int rows = n_obs;
    int cols = E_int_matrix.cols;

    if (deviceComputation) {
        //E_int_matrix.toDevice();
        //H_int_matrix.toDevice();
        double start = omp_get_wtime();
        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x);
        if (polarisation == 1) {
            computeFieldsFromMatricesKernel<<<dimGrid, dimBlock>>>(D, E_int.z, H_int.x, H_int.y, E_int_matrix.z, H_int_matrix.x, H_int_matrix.y, rows, cols);
        }
        else if (polarisation == 2) {
            computeFieldsFromMatricesKernel<<<dimGrid, dimBlock>>>(D, H_int.z, E_int.x, E_int.y, H_int_matrix.z, E_int_matrix.x, E_int_matrix.y, rows, cols);
        }
        cudaCheckError();

        

        cudaDeviceSynchronize();
        double end = omp_get_wtime();
        //printf("time = %e\n",end-start);
        //E_int.toHost();
        //H_int.toHost();

    }
    else {
        if (polarisation == 1) {
            #pragma omp parallel for
            for (int k = 0; k < rows; k++) {
                double Ez_real, Ez_imag, Hx_real, Hx_imag, Hy_real, Hy_imag;
                Ez_real = 0.0; Hx_real = 0.0; Hy_real = 0.0;
                Ez_imag = 0.0; Hx_imag = 0.0; Hy_imag = 0.0;
                
                for (int j = 0; j < cols; j++) {
                    // Get factor
                    double D_real = D.getHostRealValue(j);
                    double D_imag = D.getHostImagValue(j);

                    // Computing real values
                    Ez_real += E_int_matrix.z.getHostRealValue(k,j) * D_real;
                    Ez_real -= E_int_matrix.z.getHostImagValue(k,j) * D_imag;
                    Hx_real += H_int_matrix.x.getHostRealValue(k,j) * D_real;
                    Hx_real -= H_int_matrix.x.getHostImagValue(k,j) * D_imag;
                    Hy_real += H_int_matrix.y.getHostRealValue(k,j) * D_real;
                    Hy_real -= H_int_matrix.y.getHostImagValue(k,j) * D_imag;

                    // Computing imagninary values
                    Ez_imag += E_int_matrix.z.getHostRealValue(k,j) * D_imag;
                    Ez_imag += E_int_matrix.z.getHostImagValue(k,j) * D_real;
                    Hx_imag += H_int_matrix.x.getHostRealValue(k,j) * D_imag;
                    Hx_imag += H_int_matrix.x.getHostImagValue(k,j) * D_real;
                    Hy_imag += H_int_matrix.y.getHostRealValue(k,j) * D_imag;
                    Hy_imag += H_int_matrix.y.getHostImagValue(k,j) * D_real;
            
                }

                E_int.z.setHostRealValue(k, Ez_real);
                H_int.x.setHostRealValue(k, Hx_real);
                H_int.y.setHostRealValue(k, Hy_real);
                E_int.z.setHostImagValue(k, Ez_imag);
                H_int.x.setHostImagValue(k, Hx_imag);
                H_int.y.setHostImagValue(k, Hy_imag);
            }

        }
        else if (polarisation == 2) {
            #pragma omp parallel for
            for (int k = 0; k < rows; k++) {
                double Hz_real, Hz_imag, Ex_real, Ex_imag, Ey_real, Ey_imag, D_real, D_imag;
                Hz_real = 0.0; Ex_real = 0.0; Ey_real = 0.0;
                Hz_imag = 0.0; Ex_imag = 0.0; Ey_imag = 0.0;
                for (int j = 0; j < cols; j++) {
                    // Get factor
                    D_real = D.getHostRealValue(j);
                    D_imag = D.getHostImagValue(j);

                    // Computing real values
                    Hz_real += H_int_matrix.z.getHostRealValue(k,j) * D_real;
                    Hz_real -= H_int_matrix.z.getHostImagValue(k,j) * D_imag;
                    Ex_real += E_int_matrix.x.getHostRealValue(k,j) * D_real;
                    Ex_real -= E_int_matrix.x.getHostImagValue(k,j) * D_imag;
                    Ey_real += E_int_matrix.y.getHostRealValue(k,j) * D_real;
                    Ey_real -= E_int_matrix.y.getHostImagValue(k,j) * D_imag;

                    // Computing complex values
                    Hz_imag += H_int_matrix.z.getHostRealValue(k,j) * D_imag;
                    Hz_imag += H_int_matrix.z.getHostImagValue(k,j) * D_real;
                    Ex_imag += E_int_matrix.x.getHostRealValue(k,j) * D_imag;
                    Ex_imag += E_int_matrix.x.getHostImagValue(k,j) * D_real;
                    Ey_imag += E_int_matrix.y.getHostRealValue(k,j) * D_imag;
                    Ey_imag += E_int_matrix.y.getHostImagValue(k,j) * D_real;
            
                }

                H_int.z.setHostRealValue(k, Hz_real);
                E_int.x.setHostRealValue(k, Ex_real);
                E_int.y.setHostRealValue(k, Ey_real);
                H_int.z.setHostImagValue(k, Hz_imag);
                E_int.x.setHostImagValue(k, Ex_imag);
                E_int.y.setHostImagValue(k, Ey_imag);
            }

        }
    }
}




}