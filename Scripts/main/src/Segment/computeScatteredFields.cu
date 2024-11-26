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
void Segment::computeScatteredFields() {

    int rows = n_obs;
    int cols = E_scat_matrix.cols;

    if (deviceComputation) {
        //E_scat_matrix.toDevice();
        //H_scat_matrix.toDevice();
        double start = omp_get_wtime();
        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x);
        if (polarisation == 1) {
            computeFieldsFromMatricesKernel<<<dimGrid, dimBlock>>>(C, E_scat.z, H_scat.x, H_scat.y, E_scat_matrix.z, H_scat_matrix.x, H_scat_matrix.y, rows, cols);
        }
        else if (polarisation == 2) {
            computeFieldsFromMatricesKernel<<<dimGrid, dimBlock>>>(C, H_scat.z, E_scat.x, E_scat.y, H_scat_matrix.z, E_scat_matrix.x, E_scat_matrix.y, rows, cols);
        }
        cudaCheckError();

        

        cudaDeviceSynchronize();
        double end = omp_get_wtime();
        //printf("time = %e\n",end-start);
        //E_scat.toHost();
        //H_scat.toHost();

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
                    double C_real = C.getHostRealValue(j);
                    double C_imag = C.getHostImagValue(j);

                    // Computing real values
                    Ez_real += E_scat_matrix.z.getHostRealValue(k,j) * C_real;
                    Ez_real -= E_scat_matrix.z.getHostImagValue(k,j) * C_imag;
                    Hx_real += H_scat_matrix.x.getHostRealValue(k,j) * C_real;
                    Hx_real -= H_scat_matrix.x.getHostImagValue(k,j) * C_imag;
                    Hy_real += H_scat_matrix.y.getHostRealValue(k,j) * C_real;
                    Hy_real -= H_scat_matrix.y.getHostImagValue(k,j) * C_imag;

                    // Computing imagninary values
                    Ez_imag += E_scat_matrix.z.getHostRealValue(k,j) * C_imag;
                    Ez_imag += E_scat_matrix.z.getHostImagValue(k,j) * C_real;
                    Hx_imag += H_scat_matrix.x.getHostRealValue(k,j) * C_imag;
                    Hx_imag += H_scat_matrix.x.getHostImagValue(k,j) * C_real;
                    Hy_imag += H_scat_matrix.y.getHostRealValue(k,j) * C_imag;
                    Hy_imag += H_scat_matrix.y.getHostImagValue(k,j) * C_real;
            
                }

                E_scat.z.setHostRealValue(k, Ez_real);
                H_scat.x.setHostRealValue(k, Hx_real);
                H_scat.y.setHostRealValue(k, Hy_real);
                E_scat.z.setHostImagValue(k, Ez_imag);
                H_scat.x.setHostImagValue(k, Hx_imag);
                H_scat.y.setHostImagValue(k, Hy_imag);
            }

        }
        else if (polarisation == 2) {
            #pragma omp parallel for
            for (int k = 0; k < rows; k++) {
                double Hz_real, Hz_imag, Ex_real, Ex_imag, Ey_real, Ey_imag, C_real, C_imag;
                Hz_real = 0.0; Ex_real = 0.0; Ey_real = 0.0;
                Hz_imag = 0.0; Ex_imag = 0.0; Ey_imag = 0.0;
                for (int j = 0; j < cols; j++) {
                    // Get factor
                    C_real = C.getHostRealValue(j);
                    C_imag = C.getHostImagValue(j);

                    // Computing real values
                    Hz_real += H_scat_matrix.z.getHostRealValue(k,j) * C_real;
                    Hz_real -= H_scat_matrix.z.getHostImagValue(k,j) * C_imag;
                    Ex_real += E_scat_matrix.x.getHostRealValue(k,j) * C_real;
                    Ex_real -= E_scat_matrix.x.getHostImagValue(k,j) * C_imag;
                    Ey_real += E_scat_matrix.y.getHostRealValue(k,j) * C_real;
                    Ey_real -= E_scat_matrix.y.getHostImagValue(k,j) * C_imag;

                    // Computing complex values
                    Hz_imag += H_scat_matrix.z.getHostRealValue(k,j) * C_imag;
                    Hz_imag += H_scat_matrix.z.getHostImagValue(k,j) * C_real;
                    Ex_imag += E_scat_matrix.x.getHostRealValue(k,j) * C_imag;
                    Ex_imag += E_scat_matrix.x.getHostImagValue(k,j) * C_real;
                    Ey_imag += E_scat_matrix.y.getHostRealValue(k,j) * C_imag;
                    Ey_imag += E_scat_matrix.y.getHostImagValue(k,j) * C_real;
            
                }

                H_scat.z.setHostRealValue(k, Hz_real);
                E_scat.x.setHostRealValue(k, Ex_real);
                E_scat.y.setHostRealValue(k, Ey_real);
                H_scat.z.setHostImagValue(k, Hz_imag);
                E_scat.x.setHostImagValue(k, Ex_imag);
                E_scat.y.setHostImagValue(k, Ey_imag);
            }

        }
    }
}




}